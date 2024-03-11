import asyncio
import concurrent.futures
import re
from math import ceil
from pathlib import Path
from typing import List

import pandas as pd
from neo4j import AsyncGraphDatabase
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    SpinnerColumn,
    TextColumn,
    TimeElapsedColumn,
)
from rich.table import Table

from socialgene.compare_proteins.hmmer import CompareDomains
from socialgene.config import env_vars
from socialgene.neo4j.neo4j import GraphDriver
from socialgene.search.base import SearchBase
from socialgene.utils.logging import CONSOLE, log

progress_bar = Progress(
    TextColumn("Ingesting target clusters from database..."),
    TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
    BarColumn(),
    MofNCompleteColumn(),
    TextColumn("• Time elapsed "),
    TimeElapsedColumn(),
)


def _find_similar_proteins(
    domain_list, frac: float = 0.75, only_culture_collection: bool = False
):
    """
    The function `_find_similar_proteins` is a synchronous function that queries a Neo4j graph database to
    find similar proteins based on protein domain information.

    Args:
      protein_domain_dict (Dict[List[str]]): The `protein_domain_dict` parameter is a dictionary where
    the keys are protein uids and the values are lists of domain uids.
      frac (float): The `frac` parameter is a float value that represents the fraction of protein
    domains that need to match in order for a protein to be considered similar. By default, it is set to
    0.75, meaning that at least 75% of the protein domains need to match

    Returns:
      The function `_find_sim_protein` returns a Pandas DataFrame containing the results of the query
    executed in the Neo4j database. The DataFrame has columns `assembly_uid`, `nucleotide_uid`,
    `target`, `n_start`, and `n_end`
    """
    # TODO: move async driver to reg driver class module
    # with GraphDriver() as driver:
    with GraphDriver().driver.session() as driver:
        res = driver.run(
            f"""
                WITH $domain_list AS input_protein_domains
                MATCH (prot1:protein)<-[a1:ANNOTATES]-(h0:hmm)
                WHERE h0.uid IN input_protein_domains
                WITH input_protein_domains, prot1, count(DISTINCT(h0)) as initial_count
                WHERE initial_count > size(input_protein_domains) * $frac
                MATCH (n1:nucleotide)-[e1:ENCODES]->(prot1)
                MATCH (a1:assembly)<-[:ASSEMBLES_TO]-(n1)
                {'WHERE (a1)-[:FOUND_IN]->(:culture_collection)' if only_culture_collection else ''}
                RETURN a1.uid as assembly_uid, n1.uid as nucleotide_uid, prot1.uid as target, e1.start as n_start, e1.end as n_end
                """,
            domain_list=list(domain_list),
            frac=frac,
        ).to_df()
        return res


async def _find_similar_proteins_async(
    domain_list, frac: float = 0.75, only_culture_collection: bool = False
):
    """
    The function `_find_similar_proteins_async` is an asynchronous function that queries a Neo4j graph database to
    find similar proteins based on protein domain information.

    Args:
      protein_domain_dict (Dict[List[str]]): The `protein_domain_dict` parameter is a dictionary where
    the keys are protein uids and the values are lists of domain uids.
      frac (float): The `frac` parameter is a float value that represents the fraction of protein
    domains that need to match in order for a protein to be considered similar. By default, it is set to
    0.75, meaning that at least 75% of the protein domains need to match

    Returns:
      The function `_find_sim_protein` returns a Pandas DataFrame containing the results of the query
    executed in the Neo4j database. The DataFrame has columns `assembly_uid`, `nucleotide_uid`,
    `target`, `n_start`, and `n_end`
    """
    # TODO: move async driver to reg driver class module
    async with AsyncGraphDatabase.driver(
        env_vars["NEO4J_URI"],
        auth=(
            env_vars["NEO4J_USER"],
            env_vars["NEO4J_PASSWORD"],
        ),
    ) as driver:
        res = await driver.execute_query(
            f"""
                WITH $domain_list AS input_protein_domains
                MATCH (prot1:protein)<-[a1:ANNOTATES]-(h0:hmm)
                WHERE h0.uid IN input_protein_domains
                WITH input_protein_domains, prot1, count(DISTINCT(h0)) as initial_count
                WHERE initial_count > size(input_protein_domains) * $frac
                MATCH (n1:nucleotide)-[e1:ENCODES]->(prot1)
                MATCH (a1:assembly)<-[:ASSEMBLES_TO]-(n1)
                {'WHERE (a1)-[:FOUND_IN]->(:culture_collection)' if only_culture_collection else ''}
                RETURN a1.uid as assembly_uid, n1.uid as nucleotide_uid, prot1.uid as target, e1.start as n_start, e1.end as n_end
                """,
            domain_list=list(domain_list),
            frac=frac,
        )
        res = pd.DataFrame(res.records, columns=res.keys)
        return res


sema = asyncio.BoundedSemaphore(5)


async def _find_similar_proteins_async_multiple(
    dict_of_domain_lists,
    frac: float = 0.75,
    only_culture_collection: bool = False,
):
    # create task group
    # TODO: if webserver in future this could be used to control max time of search
    async with sema:
        async with asyncio.TaskGroup() as group:
            # create and issue tasks
            tasks = {
                k: group.create_task(
                    _find_similar_proteins_async(
                        domain_list=v,
                        frac=frac,
                        only_culture_collection=only_culture_collection,
                    )
                )
                for k, v in dict_of_domain_lists.items()
            }
        # wait for all tasks to complete...
        # report all results
        # return tasks
        return pd.concat([v.result().assign(query=k) for k, v in tasks.items()])


def _find_similar_proteins_sync_multiple(
    dict_of_domain_lists,
    frac: float = 0.75,
    only_culture_collection: bool = False,
):
    results = []
    with concurrent.futures.ThreadPoolExecutor(
        max_workers=20,
    ) as executor:
        futures = {}
        for k, v in dict_of_domain_lists.items():
            future = executor.submit(
                _find_similar_proteins,
                domain_list=v,
                frac=frac,
                only_culture_collection=only_culture_collection,
            )
            futures[future] = k
        for f in concurrent.futures.as_completed(futures):
            result = f.result()
            result["query"] = futures[f]
            results.append(result)
    return pd.concat(results)


def run_search(
    dict_of_domain_lists, frac: float = 0.75, only_culture_collection: bool = False
):
    for k, v in dict_of_domain_lists.items():
        _find_similar_proteins(
            domain_list=v,
            frac=frac,
            only_culture_collection=only_culture_collection,
        )


def run_async_search(
    dict_of_domain_lists, frac: float = 0.75, only_culture_collection: bool = False
):
    return asyncio.run(
        _find_similar_proteins_async_multiple(
            dict_of_domain_lists=dict_of_domain_lists,
            frac=frac,
            only_culture_collection=only_culture_collection,
        )
    )


def run_sync_search(
    dict_of_domain_lists, frac: float = 0.75, only_culture_collection: bool = False
):
    return _find_similar_proteins_sync_multiple(
        dict_of_domain_lists=dict_of_domain_lists,
        frac=frac,
        only_culture_collection=only_culture_collection,
    )


class SearchDomains(SearchBase, CompareDomains):
    """
    Class search for similar BGCs in a SocialGene database, using domains
    Args:
        hmm_dir (str): Path to directory containing HMM profiles (used to annotate input BGC proteins if they aren't in the database).
        use_neo4j_precalc (bool): Try to pull domains from the database or annotate all input proteins with HMMER
        modscore_cutoff (float): Minimum score cutoff for a hit to be considered significant. (modified, combined Levenshtein + Jaccard)
        **kwargs: Additional keyword arguments to pass to the parent class.
    Raises:
        ValueError: If neither sg_object nor gbk_path is provided.
    """

    def __init__(
        self,
        hmm_dir: str = None,
        input=None,
        use_neo4j_precalc: bool = True,
        modscore_cutoff: float = 0.8,
        **kwargs,
    ) -> None:
        super().__init__(input=input, modscore_cutoff=modscore_cutoff, **kwargs)
        self.hmm_dir = hmm_dir
        self.outdegree_df = pd.DataFrame
        # input sg_object or gbk_path
        self._annotate(hmm_dir=self.hmm_dir, use_neo4j_precalc=use_neo4j_precalc)
        if not self._check_for_hmm_outdegree():
            self._set_hmm_outdegree()
        self._get_outdegree_per_hmm_per_protein()

    def search(self, run_async=True, **kwargs):
        dict_of_domain_lists = (
            self.outdegree_df.groupby("protein_uid", observed=True)["hmm_uid"]
            .apply(list)
            .to_dict()
        )
        log.info(
            "Searching database for proteins with similar domain content and all of the genomes those are found in"
        )
        with Progress(
            SpinnerColumn(spinner_name="aesthetic", speed=0.2),
        ) as progress:
            task = progress.add_task("Progress...", total=2)
            progress.update(task, advance=1)
            if run_async:
                self.raw_search_results_df = run_async_search(
                    dict_of_domain_lists, **kwargs
                )
                self.raw_search_results_df["query"] = self.raw_search_results_df[
                    "query"
                ].astype("category")
                progress.update(task, advance=1)
            else:
                self.raw_search_results_df = run_sync_search(
                    dict_of_domain_lists, **kwargs
                )
                self.raw_search_results_df["query"] = self.raw_search_results_df[
                    "query"
                ].astype("category")
                progress.update(task, advance=1)
        log.info(
            f"Initial search returned {len(self.raw_search_results_df):,} proteins, found in {self.raw_search_results_df.assembly_uid.nunique():,} genomes"
        )

    def _annotate(self, hmm_dir, use_neo4j_precalc: bool = True):
        log.info("Annotating input proteins domains")
        for x in Path(hmm_dir).iterdir():
            if re.search(r"\.hmm$|hmm\.gz$", str(x)):
                self.sg_object.annotate(
                    hmm_filepath=x,
                    use_neo4j_precalc=use_neo4j_precalc,
                )

    @staticmethod
    def _check_for_hmm_outdegree() -> bool:
        """
        Checks a single HMM model to see if it has an outdegree property
        """
        with GraphDriver() as db:
            res = db.run(
                """
               MATCH (h1:hmm)
               return h1.outdegree limit 1
             """
            ).peek()
        if res.value():
            return True
        else:
            return False

    @staticmethod
    def _set_hmm_outdegree() -> None:
        """
        Writes outdegree onto all hmm nodes (not just [:ANNOTATES])
        """
        with GraphDriver() as db:
            _ = db.run(
                """
                MATCH (h1: hmm)
                SET h1.outdegree = apoc.node.degree.out(h1, "ANNOTATES")
                """
            )

    def _get_outdegree_per_hmm_per_protein(self) -> pd.DataFrame:
        """
        Get the outdegree for each domain in each protein in an input sg_object
        Args:
            sg_object (SocialGene): SocialGene object
        Returns:
            DataFrame: pd.DataFrame with columns ["protein_uid", "hmm_uid", "outdegree"]
        """
        log.info("Finding the outdegree of all domains from all input proteins")
        # Get the non-redundant set of all domains from all proteins
        doms = set()
        for v in self.sg_object.proteins.values():
            for i in v.domain_vector:
                doms.add(i)
        if doms == set():
            raise ValueError("No domains in sg_object.proteins")
        doms = list(doms)
        # Retrieve the outdegree for the non-redundant set of domains
        with GraphDriver() as db:
            db_res = db.run(
                """
               WITH $doms as doms
               MATCH (h1:hmm)
               WHERE h1.uid in doms
               RETURN h1.uid as hmm_uid, h1.outdegree as outdegree
             """,
                doms=doms,
            ).to_df()
        db_res["hmm_uid"] = db_res["hmm_uid"].astype("category")
        # Create a DF of proteins and their domains
        temp = pd.DataFrame(
            data=(
                {"protein_uid": k, "hmm_uid": i}
                for k, v in self.sg_object.proteins.items()
                for i in v.domain_vector
            ),
            dtype="category",
        )
        # return a merged df with columns ["protein_uid", "hmm_uid", "outdegree"]
        self.outdegree_df = temp.merge(db_res, how="left")
        self.outdegree_df.drop_duplicates(inplace=True, ignore_index=True)

    def _filter_max_outdegree(self, max_outdegree: int = None):
        """Filter out proteins with a higher outdegree than max_outdegree

        Args:
            max_outdegree (int): HMM model annotations with an outdegree higher than this will be dropped
        """
        # The elif != Nones are to prevent an incorrect argument from just returning the the full DF, which would be unintended
        if isinstance(max_outdegree, int) or isinstance(max_outdegree, float):
            m_start = self.outdegree_df["outdegree"].sum()
            log.info(
                f"'max_outdegree' is set to {max_outdegree:,}, will remove any domains with a higher outdegree"
            )
            self.outdegree_df = self.outdegree_df[
                self.outdegree_df["outdegree"] <= max_outdegree
            ]
            log.info(
                f"'max_outdegree' reduced the total outdegree from {m_start:,} to {self.outdegree_df['outdegree'].sum():,}"
            )
        elif max_outdegree is not None:
            raise ValueError

    def _filter_max_domains_per_protein(self, max_domains_per_protein: int = None):
        """Filter out proteins with a higher outdegree than max_outdegree

        Args:
            max_domains_per_protein (int): HMM model annotations with an outdegree higher than this will be dropped
        """
        # The elif != Nones are to prevent an incorrect argument from just returning the the full DF, which would be unintended
        if isinstance(max_domains_per_protein, int) or isinstance(
            max_domains_per_protein, float
        ):
            m_start = self.outdegree_df["outdegree"].sum()
            log.info(
                f"'max_domains_per_protein' is set to {max_domains_per_protein:,}, will remove domains from proteins from highest to lowest outdegree"
            )
            self.outdegree_df = (
                self.outdegree_df.sort_values(
                    "outdegree", ascending=True, ignore_index=True
                )
                .groupby("protein_uid", observed=True)
                .head(max_domains_per_protein)
            )
            log.info(
                f"'max_domains_per_protein' reduced the total outdegree from {m_start:,} to {self.outdegree_df['outdegree'].sum():,}"
            )
        elif max_domains_per_protein is not None:
            raise ValueError

    def _filter_max_query_proteins(self, max_query_proteins: int = None):
        """Filter out proteins with a higher outdegree than max_outdegree

        Args:
            max_query_proteins (int): HMM model annotations with an outdegree higher than this will be dropped
        """
        # The elif != Nones are to prevent an incorrect argument from just returning the the full DF, which would be unintended
        if isinstance(max_query_proteins, int) or isinstance(max_query_proteins, float):
            m_start = self.outdegree_df["outdegree"].sum()
            log.info(
                f"'max_query_proteins' is set to {max_query_proteins:,}, will limit search to {max_query_proteins} of {len(self.outdegree_df['protein_uid'].unique())} input proteins"
            )
            proteins_to_keep = (
                self.outdegree_df.sort_values(
                    "outdegree", ascending=True, ignore_index=True
                )
                .groupby("protein_uid", observed=True)
                .sum("outdegree")
                .sort_values("outdegree", ascending=True, ignore_index=False)
                .head(max_query_proteins)
                .index.to_list()
            )
            self.outdegree_df = self.outdegree_df[
                self.outdegree_df["protein_uid"].isin(proteins_to_keep)
            ]
            log.info(
                f"'max_query_proteins' reduced the total outdegree from {m_start:,} to {self.outdegree_df['outdegree'].sum():,}"
            )
        elif max_query_proteins is not None:
            raise ValueError

    def _filter_scatter(self, max_query_proteins):
        """Choose a random subset of proteins to search that are spread across the length of the input BGC."""
        log.info("Choosing query proteins that span across the input BGC")
        temp = list(self.input_bgc.features_sorted_by_midpoint)
        temp = [
            temp[int(ceil(i * len(temp) / max_query_proteins))].uid
            for i in range(max_query_proteins)
        ]
        self.outdegree_df = self.outdegree_df[
            self.outdegree_df["protein_uid"].isin(temp)
        ]
        log.info("Scattering the search to proteins to and across the input BGC")
        log.info(
            f"'max_query_proteins' is set to {max_query_proteins:,}, will limit search to {max_query_proteins} of {len(self.outdegree_df['protein_uid'].unique())} input proteins"
        )

    def prioritize_input_proteins(
        self,
        max_query_proteins: float = None,
        max_domains_per_protein: int = None,
        max_outdegree: int = None,
        scatter: bool = False,
        locus_tag_bypass_list: List[str] = None,
        protein_id_bypass_list: List[str] = None,
    ):
        """Rank input proteins by how many (:hmm)-[:ANNOTATES]->(:protein) relationships will have to be traversed

        Args:
            df (pd.DataFrame): pd.DataFrame({"external_id":[], "hmm_uid":[], "outdegree":[]})
            max_query_proteins (float): Max proteins to return. If >0 and <1, will return X% of input proteins
            max_domains_per_protein (int): Max domains to retain for each individual protein (highest outdegree dropped first)
            max_outdegree (int): HMM model annotations with an outdegree higher than this will be dropped
            scatter (bool, optional): Choose a random subset of proteins to search that are spread across the length of the input BGC. Defaults to False.
            locus_tag_bypass_list (List[str], optional): List of locus tags that will bypass filtering. This is the ID found in a GenBank file "/locus_tag=" field.  Defaults to None.
            protein_id_bypass_list (List[str], optional): Less preferred than `bypass`. List of external protein IDs that will bypass filtering. This is the ID found in a GenBank file "/protein_id=" field. Defaults to None.
        Returns:
            pd.DataFrame: pd.DataFrame({"external_id":[], "hmm_uid":[], "outdegree":[]})
        """
        log.info("Prioritizing input proteins by outdegree")
        loci_protein_ids = set()
        if locus_tag_bypass_list:
            # bypass using locus tags
            loci_protein_ids = {
                i.uid
                for i in self.input_bgc.features
                if i.locus_tag in list(locus_tag_bypass_list)
            }
        if protein_id_bypass_list:
            # bypass using an external protein ids
            loci_protein_ids.update(
                {
                    i.uid
                    for i in self.input_bgc.features
                    if i.external_id in list(protein_id_bypass_list)
                }
            )
        bypass_df = self.outdegree_df[
            self.outdegree_df["protein_uid"].isin(loci_protein_ids)
        ]
        len_start = self.outdegree_df["outdegree"].sum()
        if self.outdegree_df["outdegree"].sum() == 0:
            raise ValueError(
                "None of the domains of the input proteins have connections in the database"
            )
        if (self.outdegree_df["outdegree"] == 0).any():
            prelen = len(self.outdegree_df)
            prelen_prot = len(self.outdegree_df["protein_uid"].unique())
            self.outdegree_df = self.outdegree_df[self.outdegree_df["outdegree"] > 0]
            log.info(
                f"Removed {prelen - len(self.outdegree_df)} domains from consideration (no connections in the DB), which reemoves {prelen_prot - len(self.outdegree_df['protein_uid'].unique())} proteins from consideration"
            )

        if max_query_proteins and scatter:
            self._filter_scatter(max_query_proteins)

        self._filter_max_outdegree(max_outdegree)
        self._filter_max_domains_per_protein(max_domains_per_protein)

        if max_query_proteins and not scatter:
            self._filter_max_query_proteins(max_query_proteins)

        # Add back explicitly requested input proteins/loci

        self.outdegree_df = pd.merge(self.outdegree_df, bypass_df, how="outer")
        log.info(
            f"The total outdegree was {len_start:,}; now it's {self.outdegree_df['outdegree'].sum():,}"
        )
        self.n_searched_proteins = len(self.outdegree_df["protein_uid"].unique())

    @property
    def outdegree_table(self):
        table = Table(title="Outdegree of input protein domains")
        table.add_column("Protein", justify="left", style="cyan", no_wrap=True, ratio=1)
        table.add_column(
            "Locus/Descripton", justify="left", style="cyan", no_wrap=True, ratio=1
        )
        table.add_column(
            "Unique HMM models", justify="left", style="cyan", no_wrap=True, ratio=1
        )
        table.add_column("Mean", justify="left", style="cyan", no_wrap=True, ratio=1)
        table.add_column("Min", justify="left", style="cyan", no_wrap=True, ratio=1)
        table.add_column("25%", justify="left", style="cyan", no_wrap=True, ratio=1)
        table.add_column("50%", justify="left", style="cyan", no_wrap=True, ratio=1)
        table.add_column("75%", justify="left", style="cyan", no_wrap=True, ratio=1)
        table.add_column("Max", justify="left", style="cyan", no_wrap=True, ratio=1)
        table.add_column("Sum", justify="left", style="cyan", no_wrap=True, ratio=1)
        for i in self._outdegree_table_stats():
            table.add_row(*[str(i) for i in i])
        before_width = CONSOLE.width
        CONSOLE.width = 300
        CONSOLE.print(table)
        CONSOLE.width = before_width

    def _outdegree_table_stats(self):
        temp = (
            self.outdegree_df.groupby("protein_uid", observed=True)["outdegree"]
            .describe()
            .reset_index()
        )
        temp.drop(
            "std",
            axis=1,
            inplace=True,
        )
        summed = (
            self.outdegree_df.groupby("protein_uid", observed=True)["outdegree"]
            .sum()
            .reset_index(name="sum")
        )
        temp = pd.merge(temp, summed, on="protein_uid")
        temp = temp.sort_values("sum", ascending=True, ignore_index=True)
        temp[list(["count", "mean", "min", "25%", "50%", "75%", "max", "sum"])] = temp[
            list(["count", "mean", "min", "25%", "50%", "75%", "max", "sum"])
        ].astype(int)
        d = [
            {
                "protein_uid": i.uid,
                "desc": f"{i.external_id} | {i.locus_tag if i.locus_tag else 'no-locus-tag'} | {i.description}",
            }
            for i in self.input_bgc.features
        ]
        d = pd.DataFrame(d)
        temp = pd.merge(temp, d, on="protein_uid", how="left")
        temp = temp[
            [
                "protein_uid",
                "desc",
                "count",
                "mean",
                "min",
                "25%",
                "50%",
                "75%",
                "max",
                "sum",
            ]
        ]
        return temp.values
