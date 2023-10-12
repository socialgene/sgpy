import json
from abc import ABC, abstractmethod

import pandas as pd
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    TextColumn,
    TimeElapsedColumn,
)
from textdistance import levenshtein

from socialgene.base.socialgene import SocialGene
from socialgene.compare_proteins.hmm.hmmer import CompareDomains
from socialgene.utils.logging import log

progress_bar = Progress(
    TextColumn("Pulling target cluster data from database..."),
    BarColumn(),
    MofNCompleteColumn(),
    TextColumn("• Time elapsed "),
    TimeElapsedColumn(),
    transient=True,
)


class ProcessSearchResult:
    """Handles the processing of the DataFrame returned from searching the Neo4j Database"""

    def __init__(
        self, df: pd.DataFrame, gene_clusters_must_have_x_matches: int
    ) -> None:
        """_summary_

        Args:
            df (pd.DataFrame): pd.DataFrame(columns=['assembly_uid', 'nucleotide_uid', 'target', 'n_start', 'n_end','query'])
            gene_clusters_must_have_x_matches (int): gene_clusters_must_have_x_matches
        """
        self.df = df
        self.gene_clusters_must_have_x_matches = gene_clusters_must_have_x_matches

    def process(
        self,
        assemblies_must_have_x_matches,
        nucleotide_sequences_must_have_x_matches,
        max_gap,
    ):
        # setting columns other than "query" may actually increase the mem usage
        self.df["query"] = self.df["query"].astype("category")
        # filter assemblies, nucleotide sequences based on the number of matches to unique input BGC proteins
        self.filter_assemblies(threshold=assemblies_must_have_x_matches)
        if self.df.empty:
            raise ValueError("No hits found after filtering at the assembly level")
        self.filter_nucleotides(threshold=nucleotide_sequences_must_have_x_matches)
        if self.df.empty:
            raise ValueError("No hits found after filtering at the assembly level")
        # assign clusters of proteins that aren't interrupted by a gap greater than max_gap
        self._label_clusters(max_gap=max_gap)
        # calculate the number of unique input/query protein hits within each cluster and genome/assembly
        self._calc_intrahits()
        # self._remove_genomes_without_bgcs()
        self._remove_non_bgcs()

    def _count_unique_hits_per(self, column):
        log.info(f"Counting unique hits per {column}")
        return self.df.groupby(column)["query"].nunique().reset_index()

    def filter_assemblies(self, threshold: int):
        column = "assembly_uid"
        self._filter(column=column, threshold=threshold)

    def filter_nucleotides(self, threshold: int):
        column = "nucleotide_uid"
        self._filter(column=column, threshold=threshold)

    def _filter(self, column, threshold: int):
        if column not in ["assembly_uid", "nucleotide_uid"]:
            raise ValueError(
                f"`column` argument must be 'assembly_uid' or 'nucleotide_uid', not: {column}"
            )
        log.info(f"Filtering on {column} where unique hits >= {threshold}")
        counts = self._count_unique_hits_per(column=column)
        self.df = self.df.merge(
            counts[counts["query"] >= threshold][column],
            how="inner",
            on=column,
        )

    def _calc_intrahits(self):
        log.info("Calculating unique query hits per putative BGC")
        # count the number of unique input/query protein hits within each cluster
        self.df["cluster_unique_hits"] = self.df.groupby(
            ["assembly_uid", "nucleotide_uid", "cluster"]
        )["query"].transform("nunique")

    def _sort_genes_by_start(self):
        log.info("Sorting genes by start position")
        self.df = self.df.sort_values(
            by=["n_start"], ascending=[True], na_position="first"
        ).reset_index(inplace=False, drop=True)

    def _label_clusters(self, max_gap: int = 20000):
        """Walks through the df and for each nucleotide sequence, labels proteins within gap tolerance

        Args:
            max_gap (int, optional): Breaks a "cluster" when the diff of two proteins' starts is greater than this. Defaults to 20000.

        Returns:
            pd.DataFrame: input df + cluster column (integers)
        """
        log.info(f"Grouping protein hits if less than {str(max_gap)} bp apart")
        self._sort_genes_by_start()
        self.df["cluster"] = (
            self.df.groupby("nucleotide_uid")
            .apply(lambda group: group["n_start"] > (group["n_end"].shift() + max_gap))
            .cumsum()
            .reset_index(level=[0])[0]
        )

    @property
    def bgc_indices(self) -> pd.Index:
        """Filter the results for potential BGCs

        Returns:
             pd.Index: indices of potential BGCs
        """
        # This seems overkill but will allow dropping and re-searching in the future
        return self.df[
            self.df["cluster_unique_hits"] > self.gene_clusters_must_have_x_matches
        ].index

    @property
    def bgc_df(self) -> pd.DataFrame:
        return self.df.loc[self.bgc_indices]

    @property
    def bgc_regions(self):
        return (
            self.bgc_df.groupby(["assembly_uid", "nucleotide_uid", "cluster"])
            .agg({"n_start": "min", "n_end": "max"})
            .reset_index()
        )

    def _remove_bgcs_from_df(self):
        self.bgc_df = self.bgc_df.drop(self.bgc_indices, inplace=False)

    def _remove_genomes_without_bgcs(self):
        self.bgc_df = self.bgc_df[
            self.bgc_df["assembly_uid"].isin(
                self.bgc_df.loc[self.bgc_indices]["assembly_uid"].unique()
            )
        ]

    def _remove_non_bgcs(self):
        self.df = self.bgc_df


##################################################################################################################


class SerializeToClustermap:
    def __init__(self):
        self.uid = 0
        self.uid_dict = {}
        self.obj_df = pd.DataFrame(
            columns=["clustermap_uid", "protein_uid", "nucleotide_uid"]
        )

    def _flatten_list(x):
        return sorted([item for row in x for item in row], reverse=True)

    def _get_uid(self, obj=None):
        # increment and return uid as string
        self.uid += 1
        self.uid_dict[str(self.uid)] = obj
        return str(self.uid)

    def _build(self):
        return self._clusters() | self._create_groups_dict() | self._create_links_dict()

    # good from here
    def _clusters(self):
        log.info("Creating clustermap.js clusters")
        return {
            "clusters": [
                {
                    "uid": self._get_uid(obj=self.sg_object.assemblies[k]),
                    "name": self.sg_object.assemblies[k].name,
                    "loci": self._loci(self.sg_object.assemblies[k]),
                }
                for k in self.bgc_order
            ]
        }

    def _feature(self, feature_obj):
        feat_id = self._get_uid(obj=feature_obj)
        self.obj_df = pd.concat(
            [
                self.obj_df,
                pd.DataFrame(
                    {
                        "clustermap_uid": [feat_id],
                        "protein_uid": [feature_obj.protein_hash],
                        "nucleotide_uid": [feature_obj.parent_object.uid],
                    }
                ),
            ]
        )
        return {
            "uid": feat_id,
            "label": feature_obj.protein_id,
            "names": {
                "name": feature_obj.protein_id,
                "description": feature_obj.description,
                "hash": feature_obj.protein_hash,
            },
            "start": feature_obj.start,
            "end": feature_obj.end,
            "strand": feature_obj.strand,
        }

    def _locus(self, locus_name, locus_obj):
        return {
            "uid": self._get_uid(obj=locus_obj),
            "name": locus_name,
            "genes": [self._feature(feature_obj) for feature_obj in locus_obj.features],
            "start": min([i.start for i in locus_obj.features]),
            "end": max([i.end for i in locus_obj.features]),
        }

    def _loci(self, assembly_obj):
        return [
            self._locus(locus_name=k, locus_obj=v) for k, v in assembly_obj.loci.items()
        ]

    def _get_input_info_dict(self):
        input_info_dict = {}
        for v in self.sg_object.assemblies[self.modified_input_bgc_name].loci.values():
            for i in v.features:
                input_info_dict[i.protein_hash] = f"{i.locus_tag} {i.description}"
        return input_info_dict

    def _create_groups_dict(self):
        """Create a clustermap groups"""
        temp = pd.merge(
            self.link_df,
            self.obj_df,
            left_on=["target", "nucleotide_uid"],
            right_on=["protein_uid", "nucleotide_uid"],
            how="inner",
        )
        groups = temp.groupby(["query"])["clustermap_uid"].apply(list)
        nucleotide_uids = [
            i.uid
            for i in self.sg_object.assemblies[
                self.modified_input_bgc_name
            ].loci.values()
        ]
        groups = pd.merge(
            groups,
            self.obj_df[self.obj_df["nucleotide_uid"].isin(nucleotide_uids)],
            left_on=["query"],
            right_on=["protein_uid"],
            how="inner",
        )
        input_info_dict = self._get_input_info_dict()
        return {
            "groups": [
                {
                    "uid": self._get_uid(
                        obj=input_info_dict.get(i.protein_uid, "None")
                    ),
                    "label": input_info_dict.get(i.protein_uid, "None"),
                    "genes": [i.clustermap_uid_y] + i.clustermap_uid_x,
                }
                for x, i in groups.iterrows()
            ]
        }

    def _create_links_dict(self):
        log.info("Creating clustermap.js links")
        temp = pd.merge(
            self.link_df,
            self.obj_df,
            left_on=["target", "nucleotide_uid"],
            right_on=["protein_uid", "nucleotide_uid"],
            how="inner",
        )
        temp = temp.drop(
            labels=[
                "target",
                "protein_uid",
                "cluster",
                "assembly_uid",
                "nucleotide_uid",
            ],
            axis=1,
        )
        nucleotide_uids = [
            i.uid
            for i in self.sg_object.assemblies[
                self.modified_input_bgc_name
            ].loci.values()
        ]
        temp = pd.merge(
            temp,
            self.obj_df[self.obj_df["nucleotide_uid"].isin(nucleotide_uids)],
            left_on=["query"],
            right_on=["protein_uid"],
            how="inner",
        )
        temp = temp.drop(labels=["nucleotide_uid", "protein_uid", "query"], axis=1)
        return {
            "links": [
                {
                    "uid": self._get_uid(obj=None),
                    "target": {"uid": i.clustermap_uid_x, "name": i.clustermap_uid_x},
                    "query": {"uid": i.clustermap_uid_y, "name": i.clustermap_uid_y},
                    "identity": i.score,
                }
                for x, i in temp.iterrows()
            ]
        }

    def write_clustermap(self, outpath):
        log.info(f"Writing clustermap.js output to: {outpath}")
        with open(outpath, "w") as outfile:
            json.dump(
                self._build(),
                outfile,
            )


##################################################################################################################
class SearchBase(ABC, SerializeToClustermap):
    def __init__(
        self,
        gene_clusters_must_have_x_matches: int,
        assemblies_must_have_x_matches: int,
        nucleotide_sequences_must_have_x_matches: int,
        modscore_cutoff: float = 0.8,
    ) -> None:
        super().__init__()
        self.gene_clusters_must_have_x_matches = gene_clusters_must_have_x_matches
        self.assemblies_must_have_x_matches = assemblies_must_have_x_matches
        self.nucleotide_sequences_must_have_x_matches = (
            nucleotide_sequences_must_have_x_matches
        )
        self.input_bgc_id = None
        self.modscore_cutoff = modscore_cutoff
        self.sg_object = SocialGene()
        self.initital_search_df = pd.DataFrame()
        self.search_result = pd.DataFrame()
        self.bgc_order = pd.DataFrame()
        # links to draw between the input BGC and the target BGCs
        self.link_df = pd.DataFrame()
        # defines which input protein a target protein should be grouped with
        self.group_df = pd.DataFrame()

    @abstractmethod
    def search_db(self):
        # To be implemented in the child class
        # This function should set self.initital_search_df to a pd.DataFrame:
        # pd.DataFrame(columns=['assembly_uid', 'nucleotide_uid', 'target', 'n_start', 'n_end','query'])
        ...

    def _modify_input_bgc_name(self):
        # modify input BGC assembly name so if it also exists in the DB it won't be overridden
        self.input_bgc_id = list(self.sg_object.assemblies.keys())[0]
        self.modified_input_bgc_name = f"socialgene_query_{self.input_bgc_id}"
        self.sg_object.assemblies[
            self.modified_input_bgc_name
        ] = self.sg_object.assemblies.pop(self.input_bgc_id)
        self.sg_object.assemblies[
            self.modified_input_bgc_name
        ].uid = self.modified_input_bgc_name

    def _input_bgc_info(self):
        log.info(
            f"Input BGC has {len(self.sg_object.proteins)} proteins and/or psuedogenes"
        )

    def read_sg_object(self, sg_object: SocialGene):
        self.sg_object = sg_object
        self._modify_input_bgc_name()
        self._input_bgc_info()

    def read_input_bgc(self, gbk_path: str):
        # parse input BGC
        self.sg_object.parse(gbk_path)
        self._modify_input_bgc_name()
        self._input_bgc_info()

    def process_results(
        self,
    ):
        if self.initital_search_df.empty:
            log.warning("Stopping, no search results to process")
            return
        log.info(
            f"Starting with an initial {self.initital_search_df.assembly_uid.nunique()} putative BGCs"
        )
        self._process_search_results()
        log.info(
            f"Pulling in filtered list of {self.initital_search_df.assembly_uid.nunique()} putative BGCs"
        )
        self._fill_sg_with_hits()
        log.info(f"Filtered to {len(self.sg_object.assemblies):,} putative BGCs")
        self._create_links()
        self._choose_group()
        self._rank_order_bgcs()

    def _rank_order_bgcs(self):
        # TODO: should be based off of create_links() output not bgc_df
        """Sorts assembly IDs by comparing the levenshtein distance of ordered query proteins (forward and reverse)"""
        modified_input_bgc_name = [
            i for i in self.sg_object.assemblies if i.startswith("socialgene_query_")
        ][0]
        input_bgc_nuc_id = list(
            self.sg_object.assemblies[modified_input_bgc_name].loci.keys()
        )[0]
        input_protein_order = [
            i.protein_hash
            for i in self.sg_object.assemblies[modified_input_bgc_name]
            .loci[input_bgc_nuc_id]
            .features_sorted_by_midpoint
        ]
        temp_dict = {}
        for i in self.search_result.bgc_df["assembly_uid"].unique():
            temp_list = self.search_result.bgc_df[
                self.search_result.bgc_df["assembly_uid"] == i
            ]["query"].to_list()
            forward = levenshtein(input_protein_order, temp_list)
            temp_list.reverse()
            reverse = levenshtein(input_protein_order, temp_list)
            temp_dict[i] = min(forward, reverse)
        self.bgc_order = [self.modified_input_bgc_name]
        self.bgc_order.extend(
            [k for k, v in sorted(temp_dict.items(), key=lambda item: item[1])]
        )

    def _process_search_results(
        self,
    ):
        if isinstance(self.search_result, pd.DataFrame):
            self.search_result = ProcessSearchResult(
                df=self.initital_search_df,
                gene_clusters_must_have_x_matches=self.gene_clusters_must_have_x_matches,
            )
            self.search_result.process(
                assemblies_must_have_x_matches=self.assemblies_must_have_x_matches,
                nucleotide_sequences_must_have_x_matches=self.nucleotide_sequences_must_have_x_matches,
                max_gap=self.max_gap,
            )
        elif isinstance(self.search_result, ProcessSearchResult):
            log.debug("Search results already processed")
        else:
            raise ValueError(
                f"search_result must be a pd.DataFrame or ProcessSearchResult, not {type(self.search_result)}"
            )

    def _fill_sg_with_hits(self):
        import time

        now = time.time()
        with progress_bar as pg:
            task = pg.add_task(
                "Adding best hits...", total=len(self.search_result.bgc_regions)
            )
            for index, result in self.search_result.bgc_regions.filter(
                ["nucleotide_uid", "n_start", "n_end"]
            ).iterrows():
                self.sg_object.fill_given_locus_range(
                    locus_uid=result.iloc[0],
                    start=result.iloc[1] - 5000,
                    end=result.iloc[2] + 5000,
                )
                pg.update(task, advance=1)
        # Drop likely cross-origin proteins
        for ak, av in self.sg_object.assemblies.items():
            for nk, nv in av.loci.items():
                [i for i in nv.features if abs(i.end - i.start) < 100000]
                self.sg_object.assemblies[ak].loci[nk].features = {
                    i for i in nv.features if abs(i.end - i.start) < 100000
                }
        self.sg_object.annotate_proteins_with_neo4j(
            protein_hash_ids=None, annotate_all=True, progress=False
        )
        then = time.time()
        log.warning(f"Time to fill: {int(then - now)} seconds")

    def _prune_links(self, df: pd.DataFrame) -> pd.DataFrame:
        """Takes query to target bgc protein comparison and filters "best" hits

        Args:
            df (pd.Dataframe): pd.DataFrame({"query": [], "target": [], "score": []})

        Returns:
            pd.Dataframe: pd.DataFrame({"query": [], "target": [], "score": []})
        """
        df_to_return = pd.DataFrame()
        # Select identical protein matches first
        df = df[df["score"] >= self.modscore_cutoff]
        identical_proteins_index = df[df["query"] == df["target"]].index
        if len(identical_proteins_index) > 0:
            df_to_return = pd.concat([df_to_return, df.loc[identical_proteins_index]])
            # Remove selected proteins
            query_proteins_to_delete = set(
                df.loc[identical_proteins_index]["query"].to_list()
            )
            target_proteins_to_delete = set(
                df.loc[identical_proteins_index]["target"].to_list()
            )
            df = df[~df["query"].isin(query_proteins_to_delete)]
            df = df[~df["target"].isin(target_proteins_to_delete)]
        # For each query protein, select the target(s) with the max score
        # this can be 1 or more target proteins
        temp = df[df["score"] == df.groupby(["query"])["score"].transform("max")]
        # Remove selected proteins
        query_proteins_to_delete = set(df.loc[temp.index]["query"].to_list())
        target_proteins_to_delete = set(df.loc[temp.index]["target"].to_list())
        df_to_return = pd.concat([df_to_return, temp])
        return df_to_return

    def _create_links(self) -> pd.DataFrame:
        """
        Loop through the found putative BGCs and compare the proteins to the input BGC

        Args:
            sg_object (SocialGene): SocialGene object with one assembly id prepended with "socialgene_query_"
            bgcs (pd.DataFrame): pd.DataFrame(columns=['assembly_uid', 'nucleotide_uid', 'cluster', 'n_start', 'n_end'])

        Returns:
            pd.DataFrame: pd.DataFrame(columns=['query', 'target', 'score', 'assembly_uid', 'nucleotide_uid', 'cluster'])
        """
        # Initialize an empty DataFrame to store the links
        link_df = pd.DataFrame()
        # Get the name of the input BGC assembly
        modified_input_bgc_name = [
            i for i in self.sg_object.assemblies if i.startswith("socialgene_query_")
        ][0]
        # Get the set of protein hashes for the input BGC
        query_proteins_hash_set = self.sg_object.assemblies[
            modified_input_bgc_name
        ].protein_hash_set
        for index, row in self.search_result.bgc_regions[
            ["assembly_uid", "nucleotide_uid", "cluster"]
        ].iterrows():
            assembly_name = row.assembly_uid
            assembly = self.sg_object.assemblies[assembly_name]
            cluster = row.cluster
            nuc = row.nucleotide_uid
            # Skip the input BGC assembly
            if assembly_name == modified_input_bgc_name:
                continue
            # Get the set of protein hashes for the current assembly
            target_proteins_hash_set = assembly.protein_hash_set
            # Compare the query proteins to the target proteins using CompareDomains
            compare_domains = CompareDomains()
            compare_domains.compare_many_to_many(
                (self.sg_object.proteins[k] for k in query_proteins_hash_set),
                (self.sg_object.proteins[k] for k in target_proteins_hash_set),
            )
            # Prune the links to keep only the "best" hits
            pruned_links = self._prune_links(compare_domains.df)
            # Add the assembly_uid to the links DataFrame
            pruned_links["assembly_uid"] = assembly_name
            pruned_links["nucleotide_uid"] = nuc
            pruned_links["cluster"] = cluster
            # Add the pruned links to the links DataFrame
            link_df = pd.concat([link_df, pruned_links])
        # Convert the assembly_uid column to a categorical data type
        link_df["assembly_uid"] = link_df["assembly_uid"].astype("category")
        self.link_df = link_df

    def _choose_group(
        self,
    ):
        self.group_df = (
            self.link_df.sort_values(by="score", ascending=False)
            .groupby(
                ["target", "assembly_uid", "nucleotide_uid", "cluster"], observed=False
            )
            .head(1)
        )

    def modify_assembly_names(self):
        # Modify names of assemblies for their appearance in clustermap
        for k, v in self.sg_object.assemblies.items():
            if v.taxonomy.genus_:
                v.name = f"{v.uid}_{v.taxonomy.genus_}_{v.metadata.all_attributes['culture_collection']}"
            else:
                v.name = f"{v.uid}_{v.metadata.all_attributes['culture_collection']}"
            try:
                n_distinct_input_proteins = len(
                    self.sg_object.assemblies[
                        self.modified_input_bgc_name
                    ].protein_hash_set
                )
                unique_queries_per_assembly = (
                    self.search_result.bgc_df[
                        self.search_result.bgc_df["assembly_uid"] == k
                    ]
                    .groupby("assembly_uid")["query"]
                    .nunique()
                    .values[0]
                )
                a_perc = round(
                    unique_queries_per_assembly / len(n_distinct_input_proteins) * 100,
                    0,
                )
                n_hits = (
                    self.search_result.bgc_dfbgc_indices[
                        self.search_result.bgc_dfbgc_indices["assembly_uid"] == k
                    ]
                    .filter(["cluster", "cluster_unique_hits"])
                    .drop_duplicates()["cluster_unique_hits"]
                    .sum()
                )
                n_perc = round(n_hits / len(n_distinct_input_proteins) * 100, 0)
                v.name = f"{v.name}_{str(n_perc)}%_{str(a_perc)}%"
            except Exception as e:
                log.debug(e)
                pass
