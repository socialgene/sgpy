import asyncio
from pathlib import Path

import pandas as pd
from neo4j import AsyncGraphDatabase
from rich.console import Console
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    SpinnerColumn,
    TextColumn,
    TimeElapsedColumn,
)
from rich.table import Table

from socialgene.base.socialgene import SocialGene
from socialgene.config import env_vars
from socialgene.hmm.hmmer import HMMER
from socialgene.neo4j.neo4j import GraphDriver
from socialgene.search.base import SearchBase
from socialgene.utils.logging import log

progress_bar = Progress(
    TextColumn("Ingesting target clusters from database..."),
    TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
    BarColumn(),
    MofNCompleteColumn(),
    TextColumn("• Time elapsed "),
    TimeElapsedColumn(),
)


async def _find_similar_protein(domain_list, frac: float = 0.75):
    """
    The function `_find_sim_protein` is an asynchronous function that queries a Neo4j graph database to
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
            """
                WITH $domain_list AS input_protein_domains
                MATCH (prot1:protein)<-[a1:ANNOTATES]-(h0:hmm)
                WHERE h0.uid IN input_protein_domains
                WITH input_protein_domains, prot1, count(DISTINCT(h0)) as initial_count
                WHERE initial_count > size(input_protein_domains) * $frac
                MATCH (n1:nucleotide)-[e1:ENCODES]->(prot1)
               // WHERE (n1)-[:ASSEMBLES_TO]->(:assembly)-[:FOUND_IN]->(:culture_collection)
                MATCH (a1:assembly)<-[:ASSEMBLES_TO]-(n1)
                RETURN a1.uid as assembly_uid, n1.uid as nucleotide_uid, prot1.uid as target, e1.start as n_start, e1.end as n_end
                """,
            domain_list=list(domain_list),
            frac=frac,
        )
        res = pd.DataFrame(res.records, columns=res.keys)
        return res


async def _find_similar_protein_multiple(dict_of_domain_lists, frac: float = 0.75):
    # create task group
    # TODO: if webserver in future this could be used to control max time of search
    async with asyncio.TaskGroup() as group:
        # create and issue tasks
        tasks = {
            k: group.create_task(_find_similar_protein(domain_list=v, frac=frac))
            for k, v in dict_of_domain_lists.items()
        }
    # wait for all tasks to complete...
    # report all results
    # return tasks
    return pd.concat([v.result().assign(query=k) for k, v in tasks.items()])


def run_async_search(dict_of_domain_lists, frac: float = 0.75):
    return asyncio.run(
        _find_similar_protein_multiple(
            dict_of_domain_lists=dict_of_domain_lists, frac=frac
        )
    )


class SearchDomains(SearchBase):
    def __init__(
        self,
        hmm_dir: str = None,
        sg_object: SocialGene = None,
        gbk_path: str = None,
        max_gap: int = 20000,
        **kwargs,
    ) -> None:
        super().__init__(**kwargs)
        self.hmm_dir = hmm_dir
        self.domain_outdegree_df = pd.DataFrame
        self.gbk_path = gbk_path
        self.max_gap = max_gap
        # input sg_object or gbk_path
        if sg_object:
            self.read_sg_object(sg_object)
        elif gbk_path:
            self.read_input_bgc(self.gbk_path)
        else:
            raise ValueError("Must provide either sg_object or gbk_path")
        self._annotate(hmm_dir=self.hmm_dir)
        if not self._check_for_hmm_outdegree():
            self._set_hmm_outdegree()
        self._get_outdegree_per_hmm_per_protein()

    def search_db(self):
        dict_of_domain_lists = (
            self.domain_outdegree_df.groupby("protein_uid", observed=True)["hmm_uid"]
            .apply(list)
            .to_dict()
        )
        log.info("Searching database for proteins with similar domain content")
        with Progress(
            SpinnerColumn(spinner_name="aesthetic", speed=0.2),
        ) as progress:
            task = progress.add_task("Progress...", total=2)
            progress.update(task, advance=1)
            self.initital_search_df = run_async_search(dict_of_domain_lists)
            progress.update(task, advance=1)
        log.info(
            f"Initial search returned {len(self.initital_search_df):,} proteins, found in {self.initital_search_df.assembly_uid.nunique():,} genomes"
        )

    def _annotate(self, hmm_dir):
        log.info("Annotating input proteins domains")
        if hmm_dir:
            # hmm file with cutoffs
            h1 = Path(hmm_dir, "socialgene_nr_hmms_file_with_cutoffs_1_of_1.hmm")
            # hmm file without cutoffs
            h2 = Path(hmm_dir, "socialgene_nr_hmms_file_without_cutoffs_1_of_1.hmm")
            hmm_file_list = [h1, h2]
            for i in hmm_file_list:
                HMMER().hmmpress(i)
        self.sg_object.annotate(
            hmm_directory=hmm_dir,
            use_neo4j_precalc=True,
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
        self.domain_outdegree_df = temp.merge(db_res, how="left")
        self.domain_outdegree_df.drop_duplicates(inplace=True, ignore_index=True)

    def prioritize_input_proteins(
        self,
        max_query_proteins: int = None,
        max_domains_per_protein: int = None,
        max_outdegree: int = None,
    ):
        """Rank input proteins by how many (:hmm)-[:ANNOTATES]->(:protein) relationships will have to be traversed

        Args:
            df (pd.DataFrame): pd.DataFrame({"protein_id":[], "hmm_uid":[], "outdegree":[]})
            max_query_proteins (int): Max proteins to return
            max_domains_per_protein (int): Max domains to retain for each individual protein (highest outdegree dropped first)
            max_outdegree (int): HMM model annotations with an outdegree higher than this will be dropped

        Returns:
            pd.DataFrame: pd.DataFrame({"protein_id":[], "hmm_uid":[], "outdegree":[]})
        """
        log.info("Prioritizing input proteins by outdegree")
        len_start = self.domain_outdegree_df["outdegree"].sum()
        # The elif != Nones are to prevent an incorrect argument from just returning the the full DF, which would be unintended (e.g. someone uses a float)

        if isinstance(max_outdegree, int) or isinstance(max_outdegree, float):
            self.domain_outdegree_df = self.domain_outdegree_df[
                self.domain_outdegree_df["outdegree"] <= max_outdegree
            ]
        elif max_outdegree is not None:
            raise ValueError

        if isinstance(max_domains_per_protein, int) or isinstance(
            max_domains_per_protein, float
        ):
            self.domain_outdegree_df = (
                self.domain_outdegree_df.sort_values(
                    "outdegree", ascending=True, ignore_index=True
                )
                .groupby("protein_uid", observed=True)
                .head(max_domains_per_protein)
            )
        elif max_domains_per_protein is not None:
            raise ValueError
        if isinstance(max_query_proteins, int) or isinstance(max_query_proteins, float):
            proteins_to_keep = (
                self.domain_outdegree_df.sort_values(
                    "outdegree", ascending=True, ignore_index=True
                )
                .groupby("protein_uid", observed=True)
                .sum("outdegree")
                .sort_values("outdegree", ascending=True, ignore_index=False)
                .head(max_query_proteins)
                .index.to_list()
            )
            self.domain_outdegree_df = self.domain_outdegree_df[
                self.domain_outdegree_df["protein_uid"].isin(proteins_to_keep)
            ]
        elif max_query_proteins is not None:
            raise ValueError
        log.info(
            f"domain_outdegree_df was {len_start:,}, now: {self.domain_outdegree_df['outdegree'].sum():,}"
        )

    def outdegree_table(self):
        table = Table(title="Outdegree of input proteins")
        table.add_column("Protein", justify="left", style="cyan", no_wrap=True, ratio=1)
        table.add_column(
            "Distinct HMMs", justify="left", style="cyan", no_wrap=True, ratio=1
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
        console = Console()
        console.print(table)

    def _outdegree_table_stats(self):
        temp = (
            self.domain_outdegree_df.groupby("protein_uid", observed=True)["outdegree"]
            .describe()
            .reset_index()
        )
        temp.drop(
            "std",
            axis=1,
            inplace=True,
        )
        summed = (
            self.domain_outdegree_df.groupby("protein_uid", observed=True)["outdegree"]
            .sum()
            .reset_index(name="sum")
        )
        temp = pd.merge(temp, summed, on="protein_uid")
        temp = temp.sort_values("sum", ascending=True, ignore_index=True)
        temp[list(["count", "mean", "min", "25%", "50%", "75%", "max", "sum"])] = temp[
            list(["count", "mean", "min", "25%", "50%", "75%", "max", "sum"])
        ].astype(int)
        return temp.values
