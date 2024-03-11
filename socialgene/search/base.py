import time
from abc import ABC, abstractmethod
from math import ceil
from pathlib import Path

import pandas as pd
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    TextColumn,
    TimeElapsedColumn,
)
from rich.table import Table
from textdistance import levenshtein

from socialgene.base.socialgene import SocialGene
from socialgene.clustermap.serialize import SerializeToClustermap
from socialgene.compare_gene_clusters.compare_gene_clusters import BGCComparison
from socialgene.neo4j.neo4j import GraphDriver
from socialgene.utils.logging import CONSOLE, log

progress_bar = Progress(
    TextColumn("Pulling target cluster data from database..."),
    BarColumn(),
    MofNCompleteColumn(),
    TextColumn("• Time elapsed "),
    TimeElapsedColumn(),
    transient=True,
)


class BgcComp2:
    def __init__(self, jaccard, levenshtein, modscore):
        self.jaccard = jaccard
        self.levenshtein = levenshtein
        self.modscore = modscore


def truncate_string(str_input, max_length):
    str_end = ".."
    length = len(str_input)
    if length > max_length:
        return str_input[: max_length - len(str_end)] + str_end
    else:
        return str_input


##################################################################################################################
class SearchBase(ABC):
    """
    Base class for searching for a biosynthetic gene clusters (BGCs) in a SocialGene database

    Args:
        gene_clusters_must_have_x_matches (float): Minimum number of distinct query protein matches required for a BGC to be considered a hit. <1 is a fraction of the number of query proteins, >1 is the number of query proteins.
        assemblies_must_have_x_matches (float): Minimum number of distinct query protein matches for a genome to be considered. <1 is a fraction of the number of query proteins, >1 is the number of query proteins.
        nucleotide_sequences_must_have_x_matches (float): Minimum number of distinct query protein matches required for a nucleotide sequence to be considered. <1 is a fraction of the number of query proteins, >1 is the number of query proteins.
        modscore_cutoff (float, optional): Minimum score required for a BGC to be considered a hit. Defaults to 0.8.
        break_bgc_on_gap_of(int, optional): Putative BGCs are split when matched proteins are separated by n base pairs. Defaults to 20000.
        target_bgc_padding (int, optional): Number of base pairs to add to the start and end of a target BGC. Defaults to 10000.


    Attributes:
        gene_clusters_must_have_x_matches (float): Minimum number of distinct query protein matches required for a BGC to be considered a hit. <1 is a fraction of the number of query proteins, >1 is the number of query proteins.
        assemblies_must_have_x_matches (float): Minimum number of distinct query protein matches for a genome to be considered. <1 is a fraction of the number of query proteins, >1 is the number of query proteins.
        nucleotide_sequences_must_have_x_matches (float): Minimum number of distinct query protein matches required for a nucleotide sequence to be considered. <1 is a fraction of the number of query proteins, >1 is the number of query proteins.
        input_bgc_id (str): ID of the input biosynthetic gene cluster (BGC).
        input_assembly (Assembly): Assembly object of the input BGC.
        modscore_cutoff (float): Minimum score required for a BGC to be considered a hit.
        break_bgc_on_gap_of (int): Number of base pairs to break a BGC on if a gap is greater than this value.
        sg_object (SocialGene): SocialGene object containing the database.
        raw_search_results_df (pd.DataFrame): Initial search results as a pandas DataFrame.
        sorted_bgcs (pd.DataFrame): Order of the BGCs in the search results as a pandas DataFrame.
        link_df (pd.DataFrame): Linkage information between BGCs in the search results as a pandas DataFrame.
        group_df (pd.DataFrame): Grouping information for BGCs in the search results as a pandas DataFrame.
        n_searched_proteins (int): Number of proteins searched in the database.
        target_bgc_padding (int): Number of base pairs to add to the start and end of a target BGC.
        working_search_results_df (pd.DataFrame): Working search results as a pandas DataFrame.

    Methods:
        input (): SocialGene object containing a single BGC to search (optional, must provide this or a gbk_path); or Path to GenBank file containing a single BGC to search (optional, must provide this or a sg_object).
        search(): Abstract method to be implemented in child classes. Searches the database for BGCs.
        _modify_input_bgc_name(): Modifies the name of the input BGC to avoid overriding existing BGCs in the database.
        _input_bgc_info(): Logs information about the input BGC.
        read_sg_object(sg_object: SocialGene): Reads a SocialGene object as the database.
        read_input_bgc(gbk_path: str): Reads an input BGC from a GenBank file.
        process_results(): Processes the search results to filter and group BGCs.
    """

    def __init__(
        self,
        input,
        gene_clusters_must_have_x_matches: float,
        assemblies_must_have_x_matches: float,
        nucleotide_sequences_must_have_x_matches: float,
        modscore_cutoff: float = 0.8,
        target_bgc_padding=10000,
        break_bgc_on_gap_of: int = 20000,
    ) -> None:
        super().__init__()
        self.gene_clusters_must_have_x_matches = gene_clusters_must_have_x_matches
        self.assemblies_must_have_x_matches = assemblies_must_have_x_matches
        self.nucleotide_sequences_must_have_x_matches = (
            nucleotide_sequences_must_have_x_matches
        )
        self.input_bgc_id = None
        self.input_assembly = None
        self.modscore_cutoff = modscore_cutoff
        self.break_bgc_on_gap_of = break_bgc_on_gap_of
        self.sg_object = SocialGene()
        self.raw_search_results_df = pd.DataFrame()
        self.sorted_bgcs = None
        self.link_df = pd.DataFrame()
        self.group_df = pd.DataFrame()
        self.n_searched_proteins = None
        self.target_bgc_padding = target_bgc_padding
        self.working_search_results_df = pd.DataFrame()
        self.outdegree_df = pd.DataFrame()
        self.bgc_df = pd.DataFrame()
        self.primary_bgc_regions = pd.DataFrame()

        if (
            isinstance(input, SocialGene)
            or isinstance(input, str)
            or isinstance(input, Path)
        ):
            self.read_input_bgc(input)
            self.gbk_path = input
        else:
            raise ValueError("Must provide either sg_object or gbk_path")
        self.n_searched_proteins = len(self.sg_object.proteins)
        self._compare_two_gene_clusters_score = BgcComp2

    @abstractmethod
    def prioritize_input_proteins(
        self,
        scatter: bool = False,
    ):
        # To be implemented in the child class
        # scatter (bool, optional): Choose a random subset of proteins to search that are spread across the length of the input BGC. Defaults to False.
        ...

    @abstractmethod
    def search(self):
        # To be implemented in the child class
        # This function should set self.raw_search_results_df to a pd.DataFrame:
        # pd.DataFrame(columns=['assembly_uid', 'nucleotide_uid', 'target', 'n_start', 'n_end','query'])
        ...

    @property
    @abstractmethod
    def outdegree_table(self):
        # prints info about the prioritized proteins for the input BGC
        ...

    @abstractmethod
    def _outdegree_table_stats(self):
        # generates table for outdegree_table()
        ...

    def _modify_input_bgc_name(self):
        # modify input BGC assembly name so if it also exists in the DB it won't be overridden
        self.input_bgc_id = list(self.sg_object.assemblies.keys())[0]
        self.modified_input_bgc_name = f"socialgene_query_{self.input_bgc_id}"
        self.sg_object.assemblies[self.modified_input_bgc_name] = (
            self.sg_object.assemblies.pop(self.input_bgc_id)
        )
        self.sg_object.assemblies[self.modified_input_bgc_name].uid = (
            self.modified_input_bgc_name
        )
        self.input_assembly = self.sg_object.assemblies[self.modified_input_bgc_name]

    def create_input_bgcs(self):
        for locus in self.input_assembly.loci.values():
            locus.add_bgcs_by_feature(features=locus.features)
        self.input_bgc = list(self.input_assembly.gene_clusters)[0]

    def _input_bgc_info(self):
        log.info(
            f"Input BGC has {len(self.sg_object.proteins)} proteins and/or pseudogenes"
        )

    def read_input_bgc(self, input: str):
        # parse input BGC
        if isinstance(input, SocialGene):
            self.sg_object = input
        else:
            self.sg_object.parse(input)
        self._modify_input_bgc_name()
        self._input_bgc_info()
        self.create_input_bgcs()

    def filter(self, drop_raw_search_results_df=False):
        if self.raw_search_results_df.empty:
            log.warning("Stopping, no search results to process")
            return
        log.info(
            f"Starting with matches across {self.raw_search_results_df.assembly_uid.nunique()} genomes"
        )
        self.working_search_results_df = self.raw_search_results_df
        if drop_raw_search_results_df:
            del self.raw_search_results_df
        # filter assemblies, nucleotide sequences based on the number of matches to unique input BGC proteins
        self._filter_assemblies(threshold=self.assemblies_must_have_x_matches)
        if self.working_search_results_df.empty:
            ValueError("No hits found after filtering at the assembly level")
        self._filter_nucleotides(
            threshold=self.nucleotide_sequences_must_have_x_matches
        )
        if self.working_search_results_df.empty:
            raise ValueError("No hits found after filtering at the assembly level")

    # def cluster(self):
    #     # assign clusters of proteins that aren't interrupted by a gap greater than break_bgc_on_gap_of
    #     self._label_clusters()

    def _primary_bgc_regions(self, limiter=None):
        return self._collapse_cluster(
            self._filter_clusters(
                df=self.working_search_results_df,
                threshold=self.gene_clusters_must_have_x_matches,
                limiter=limiter,
            )
        )

    def annotate(self):
        try:
            self._create_links()
            self._choose_group()
            self._rank_order_bgcs()
        except Exception as e:
            log.warning(e)

    def _compare_two_gene_clusters(self, group, q_len):
        # levenshtein similarity of query proteins to target proteins
        # compares the order of query proteins to query proteins ordered by rbh to target
        only_matched = group.dropna(subset=["query"])
        a = levenshtein(
            list(only_matched.sort_values(["t_start"], ascending=True)["query"]),
            list(only_matched.sort_values(["q_start"], ascending=True)["query"]),
        ) / len(only_matched)
        # levenshtein similarity of query proteins to target proteins in reverse order
        b = levenshtein(
            list(only_matched.sort_values(["t_start"], ascending=True)["query"]),
            list(only_matched.sort_values(["q_start"], ascending=False)["query"]),
        ) / len(only_matched)
        lev = a if a < b else b
        lev = 1 - lev
        jac = (len(group[group["target_feature"].notnull()])) / len(
            group["query_feature"].unique()
        )

        return self._compare_two_gene_clusters_score(
            jac, lev, (jac * 2) + lev - (abs(len(group) - q_len) / q_len)
        )

        # return self._compare_two_gene_clusters_score(
        #     jac, lev, (jac * 2) + lev - (abs(len(group) - q_len) / q_len)
        # )

    def _compare_all_gcs_to_input_gc(self):
        # create a 3-column dataframe that represents the input gene cluster's features (proteins)
        q_obj = pd.DataFrame(
            [
                {
                    "query_feature": i,
                    "q_start": i.start,
                    "strand": i.strand,
                }
                for i in list(self.input_bgc.features_sorted_by_midpoint)
            ]
        )
        for k, group in self.link_df.groupby(
            ["query_gene_cluster", "target_gene_cluster"], sort=False
        ):
            score = self._compare_two_gene_clusters(
                pd.merge(
                    q_obj,
                    group,
                    on="query_feature",
                    how="outer",
                    suffixes=("", "_delme"),
                ),
                len(self.input_bgc.features),
            )
            yield {
                "query_gene_cluster": k[1],
                "target_gene_cluster": k[0],
                "jaccard": score.jaccard,
                "levenshtein": score.levenshtein,
                "modscore": score.modscore,
            }

    def _compare_bgcs_by_jaccard_and_levenshtein(
        self,
    ):
        """Sorts assembly IDs by comparing the levenshtein distance of ordered query proteins (forward and reverse)"""
        return pd.DataFrame(self._compare_all_gcs_to_input_gc()).sort_values(
            by="modscore", ascending=False
        )

    def _compare_bgcs_by_median_bitscore(
        self,
    ):
        temp = self.link_df[["query_gene_cluster", "target_gene_cluster", "score"]]
        return (
            temp.groupby(["query_gene_cluster", "target_gene_cluster"])["score"]
            .median()
            .sort_values(ascending=False)
            .reset_index()
        )

    def _fraction_matched(self):
        df = (
            self.link_df.groupby(["target_gene_cluster"])["query"]
            .nunique()
            .apply(lambda x: x / len(self.input_bgc.proteins))
            .reset_index()
        )
        df.columns = ["gene_cluster", "fraction_matched"]
        return df

    def _rank_order_bgcs(self, threshold):
        """Sorts assemblies by comparing the jaccard and levenshtein distance of ordered
            query proteins (forward and reverse) and breaking ties with the median bitscore of RBH between BGCs

        Args:
            threshold (float): The minimum jaccard score required for a BGC to be considered a hit. Think of this as the number of input bgc proteins a putative bgc must have to be considered a hit. (0 < threshold < 1)
        Returns:
            list[Assembly]: LIst of SocialGene Assembly objects
        """
        temp = pd.merge(
            self._compare_bgcs_by_jaccard_and_levenshtein(),
            self._compare_bgcs_by_median_bitscore(),
            left_on="query_gene_cluster",
            right_on="target_gene_cluster",
            how="inner",
        )
        temp.drop(
            ["query_gene_cluster_y", "target_gene_cluster_y"], axis=1, inplace=True
        )
        temp = temp.sort_values(by=["modscore", "score"], ascending=False)
        secondary = temp[temp["jaccard"] < threshold]
        primary = temp[temp["jaccard"] >= threshold]
        for i, row in secondary.iterrows():
            # remove the gene cluster from the sg_object
            z = row.query_gene_cluster_x
            z.parent.gene_clusters = [i for i in z.parent.gene_clusters if i != z]
        primary["query_gene_cluster_x_assembly"] = primary[
            "query_gene_cluster_x"
        ].apply(lambda x: x.parent.parent.uid)
        primary.drop_duplicates(subset=["query_gene_cluster_x_assembly"], inplace=True)
        return primary["query_gene_cluster_x"].to_list()

    def _bgc_regions_to_sg_object(self, df):
        now = time.time()
        log.info(
            f"Pulling data from the database for {df.nucleotide_uid.nunique()} putative BGCs"
        )
        with progress_bar as pg:
            task = pg.add_task("Adding best hits...", total=len(df))
            for index, df_row in df.iterrows():
                # pull in info from database to sg_object
                _ = self.sg_object.fill_given_locus_range(
                    locus_uid=df_row.loc["nucleotide_uid"],
                    start=df_row.loc["n_start"] - self.target_bgc_padding,
                    end=df_row.loc["n_end"] + self.target_bgc_padding,
                )
                # add bgc to locus
                pg.update(task, advance=1)
        self.sg_object._drop_all_cross_origin()

        then = time.time()
        log.info(f"Time to fill: {int(then - now)} seconds")

    def _create_link_df(
        self, query_gene_cluster, target_gene_cluster, tool="blastp", **kwargs
    ):
        """read the args, things are backward b/c context is searching the found bgc against the input bgc

        Args:
            query_gene_cluster (_type_): the database putative gene cluster
            target_gene_cluster (_type_): the user input gene cluster

        Returns:
            _type_: _description_
        """

        comparator = BGCComparison(
            tool=tool,
        )
        protein_comparisons_df = comparator.compare(
            query_gene_cluster.protein_iter, target_gene_cluster.protein_iter, **kwargs
        )
        # Turn the protein uid in the blastp etc output into protein objects from the cluster sg_objects
        protein_comparisons_df["query_protein"] = protein_comparisons_df["query"].apply(
            lambda x: query_gene_cluster.proteins[x]
        )
        protein_comparisons_df["target_protein"] = protein_comparisons_df[
            "target"
        ].apply(lambda x: target_gene_cluster.proteins[x])
        # expand to get one row per feature per protein (one non-redunant protein can have multiple features)
        q_obj = pd.DataFrame(
            [
                {
                    "query_gene_cluster": query_gene_cluster,
                    "query_feature": i,
                    "query_protein": i.protein,
                    "q_start": i.start,
                    "strand": i.strand,
                }
                for i in list(query_gene_cluster.features_sorted_by_midpoint)
            ]
        )
        # expand to get one row per feature per protein (one non-redunant protein can have multiple features)
        t_obj = pd.DataFrame(
            [
                {
                    "target_gene_cluster": target_gene_cluster,
                    "target_feature": i,
                    "target_protein": i.protein,
                    "t_start": i.start,
                    "strand": i.strand,
                }
                for i in list(target_gene_cluster.features_sorted_by_midpoint)
            ]
        )
        q_temp = pd.merge(protein_comparisons_df, q_obj, on="query_protein", how="left")
        return (
            pd.merge(q_temp, t_obj, on="target_protein", how="left")
            .sort_values("t_start")
            .rename(columns={comparator.tool.score_column: "score"}, inplace=False)
        )

    def _create_links(self, tool="hmmer", **kwargs) -> pd.DataFrame:
        """
        Loop through gene_cluster compare the proteins to the input BGC

        Returns:
        pd.DataFrame: A DataFrame containing the links between the input BGC and the putative BGCs.
        The DataFrame has the following columns:
            - query: The protein hash of the input BGC protein.
            - target: The protein hash of the putative BGC protein.
            - score: The score of the link between the input BGC protein and the putative BGC protein.
            - assembly_uid: The assembly UID of the putative BGC.
            - nucleotide_uid: The nucleotide UID of the putative BGC.
            - cluster: The cluster of the putative BGC.
        """
        # Initialize an empty DataFrame to store the links
        log.info("Start: Creating links")
        link_df = pd.DataFrame()
        progress_bar = Progress(
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            MofNCompleteColumn(),
            TextColumn("• Time elapsed "),
            TimeElapsedColumn(),
            transient=True,
        )
        log.info(
            f"Finding reciprocal best hits; protein similarity via {BGCComparison(tool=tool).tool.name}"
        )
        with progress_bar as pg:
            task = pg.add_task(
                description="Creating links",
                total=len(list(self.sg_object.get_all_gene_clusters())),
            )
            for target_gene_cluster in self.sg_object.get_all_gene_clusters():
                if target_gene_cluster != self.input_bgc:
                    pg.update(
                        task,
                        description=(
                            f"\nComparing {self.input_bgc_id} to {target_gene_cluster.parent.external_id}:{min([i.start for i in target_gene_cluster.features])}-{max([i.start for i in target_gene_cluster.features])}"
                        ),
                    )
                    link_df = pd.concat(
                        [
                            link_df,
                            self._create_link_df(
                                query_gene_cluster=self.input_bgc,
                                target_gene_cluster=target_gene_cluster,
                                tool=tool,
                                **kwargs,
                            ),
                        ]
                    )

                pg.update(task, advance=1, description="[blue]Complete Task")
        self.link_df = link_df
        log.info(f"Finish: Creating links; {len(link_df)} links produced")
        if self.link_df.empty:
            log.warning("No links produced")

    def _choose_group(
        self,
    ):
        log.info("Start: Assigning target BGC proteins to input BGC groups")
        if self.link_df.empty:
            log.info("Finish: Assigning target BGC proteins to input BGC groups")
            log.warning("No links to group by, no groups produced")
            return None
        self.group_df = (
            self.link_df.sort_values(by="score", ascending=False)
            .groupby(["query"], observed=False)
            .head(1)
        )
        log.info("Finish: Assigning target BGC proteins to input BGC groups")

    def add_info_to_cluster_names_for_clustermapjs(self):
        # Modify names of assemblies for their appearance in clustermap
        for k, v in self.sg_object.assemblies.items():
            if v.taxonomy.genus_:
                v.name = f"{v.uid}_{v.taxonomy.genus_}_{v.metadata.all_attributes()['culture_collection']}"
            else:
                v.name = f"{v.uid}_{v.metadata.all_attributes()['culture_collection']}"
            try:
                n_distinct_input_proteins = len(self.input_assembly.feature_uid_set)
                unique_queries_per_assembly = (
                    self.bgc_df[self.bgc_df["assembly_uid"] == k]
                    .groupby("assembly_uid")["query"]
                    .nunique()
                    .values[0]
                )
                a_perc = round(
                    unique_queries_per_assembly / len(n_distinct_input_proteins) * 100,
                    0,
                )
                n_hits = (
                    self.bgc_dfbgc_indices[self.bgc_dfbgc_indices["assembly_uid"] == k]
                    .filter(["cluster", "cluster_unique_hits"])
                    .drop_duplicates()["cluster_unique_hits"]
                    .sum()
                )
                n_perc = round(n_hits / len(n_distinct_input_proteins) * 100, 0)
                v.name = f"{v.name}_{str(n_perc)}%_{str(a_perc)}%"
            except Exception as e:
                log.debug(e)
                pass

    def _count_unique_hits_per(self, column):
        log.info(f"Counting unique hits per {column}")
        return (
            self.working_search_results_df.groupby(column)["query"]
            .nunique()
            .reset_index()
        )

    def _filter_assemblies(self, threshold: float):
        column = "assembly_uid"
        self._filter(column=column, threshold=threshold)

    def _filter_nucleotides(self, threshold: float):
        column = "nucleotide_uid"
        self._filter(column=column, threshold=threshold)

    def _filter(self, column, threshold: float):
        log.info(f"Filtering on {column} where unique hits >= {threshold}")
        if threshold > 0 and threshold < 1:
            threshold = ceil(self.n_searched_proteins * threshold)
        counts = self._count_unique_hits_per(column=column)
        self.working_search_results_df = self.working_search_results_df.merge(
            counts[counts["query"] >= threshold][column],
            how="inner",
            on=column,
        )
        log.info(
            f"{self.working_search_results_df['assembly_uid'].nunique()} assemblies, {self.working_search_results_df['nucleotide_uid'].nunique()} nucleotide sequences had {column}s with >= {threshold} unique query hits"
        )

    def _count_unique_hits_per_cluster(self, df):
        return df.groupby(["cluster"])["query"].nunique().sort_values(ascending=False)

    def _sort_genes_by_start(self):
        log.info("Sorting genes by start position")
        self.working_search_results_df = self.working_search_results_df.sort_values(
            by=["n_start"], ascending=[True], na_position="first"
        ).reset_index(inplace=False, drop=True)

    def label_clusters(
        self,
    ):
        """Walks through the df and for each nucleotide sequence, assigns genes/protein to a cluster, breaking on gap tolerance

        Args:
            break_bgc_on_gap_of (int, optional): Breaks a "cluster" when the diff of two proteins' starts is greater than this. Defaults to 20000.

        Returns:
            pd.DataFrame: input df + cluster column (integer)
        """
        log.info(
            f"Grouping protein hits if less than {str(self.break_bgc_on_gap_of)} bp apart"
        )
        self._sort_genes_by_start()
        # have to check if groupby output is a df or series and adjust accordingly
        # https://stackoverflow.com/questions/37715246/pandas-groupby-apply-behavior-returning-a-series-inconsistent-output-type
        n_groups = self.working_search_results_df["nucleotide_uid"].nunique()
        if n_groups == 1:
            self.working_search_results_df["cluster"] = (
                self.working_search_results_df["n_start"]
                - self.working_search_results_df["n_end"].shift().fillna(0).astype(int)
                > self.break_bgc_on_gap_of
            ).cumsum()
        elif n_groups > 1:
            self.working_search_results_df["cluster"] = (
                self.working_search_results_df.groupby("nucleotide_uid")
                .apply(
                    lambda group: group["n_start"]
                    > (group["n_end"].shift() + self.break_bgc_on_gap_of)
                )
                .cumsum()
                .reset_index(level=[0])[0]
            )
        else:
            raise ValueError("Huh, expected one or more groups")

    def _collapse_cluster(self, df):
        """Collapse a cluster to a single row

        Args:
            df (pd.DataFrame): input df

        Returns:
            pd.DataFrame: collapsed df
        """
        return (
            df.groupby(["assembly_uid", "nucleotide_uid", "cluster"])
            .agg({"n_start": "min", "n_end": "max"})
            .reset_index()
        )

    def _filter_clusters(self, df, threshold, limiter=None) -> pd.Index:
        """Filter the results for potential BGCs
        Args:
        Returns:
            pd.Index: indices of potential BGCs
        """
        if threshold > 0 and threshold < 1:
            threshold = ceil(self.n_searched_proteins * threshold)

        tempdf = self._count_unique_hits_per_cluster(df)
        if limiter:
            tempdf = tempdf[0:limiter]
        tempdf = tempdf[tempdf >= threshold]
        df = df.filter(items=tempdf.index, axis=0)
        return df

    def write_clustermap_json(self, outpath):
        cmap = SerializeToClustermap(
            sorted_bgcs=self.sorted_bgcs, sg_object=self.sg_object, link_df=self.link_df
        )
        cmap.write(outpath=outpath)

    def user_friendly_hit_df(self, truncate_description_to_n_chars=20):
        with GraphDriver() as db:
            db_res = db.run(
                """
               WITH $input as inputs
               unwind inputs as input
               MATCH (n1:nucleotide {uid:input.nucleotide_uid})-[e1:ENCODES {start:input.n_start}]->(p1:protein {uid:input.target})
               RETURN input.assembly_uid as t_assembly_uid,
                      n1.external_id as t_external_id,
                      e1.locus_tag as t_locus_tag,
                      e1.protein_id as t_protein_id,
                      e1.start as t_start,
                      e1.end as t_end,
                      e1.description as t_description,
                      input.query as query,
                      input.cluster as cluster
             """,
                input=self.working_search_results_df.to_dict("records"),
            ).to_df()
        query_df = pd.DataFrame(
            [
                {
                    "query": i.uid,
                    "q_description": i.description,
                    "q_external_id": i.external_id,
                }
                for i in self.input_bgc.features
            ]
        )
        db_res = pd.merge(db_res, query_df, on="query", how="left")
        db_res = db_res.drop(["query"], axis=1)
        db_res = db_res[
            [
                "t_assembly_uid",
                "t_external_id",
                "t_locus_tag",
                "cluster",
                "q_external_id",
                "t_protein_id",
                "t_start",
                "t_end",
                "t_description",
                "q_description",
            ]
        ]
        db_res = db_res.sort_values(["t_assembly_uid", "t_external_id", "t_start"])
        db_res["t_description"] = db_res["t_description"].apply(
            lambda x: truncate_string(x, truncate_description_to_n_chars)
        )
        db_res["q_description"] = db_res["q_description"].apply(
            lambda x: truncate_string(x, truncate_description_to_n_chars)
        )
        return db_res

    def rich_table(
        self,
        df,
        title="",
        n=100,
        justify="center",
        style="cyan",
        no_wrap=True,
        **kwargs,
    ):
        table = Table(title=title)
        for i in df.columns:
            table.add_column(
                i, justify=justify, style=style, no_wrap=no_wrap, ratio=1, **kwargs
            )
        ii = 0
        for i in df.values:
            table.add_row(*[str(i) for i in i])
            ii += 1
            if ii == n:
                break
        console = CONSOLE
        console.print(table)
