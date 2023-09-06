import itertools
from multiprocessing import Pool
from socialgene.clustermap.clustermap import Clustermap
from socialgene.hmm.hmmer import HMMER
from socialgene.base.socialgene import SocialGene
from pathlib import Path
from socialgene.scoring.scoring import mod_score

from socialgene.scoring.search import (
    check_for_hmm_outdegree,
    get_lowest_outdegree_model_per_protein,
    search_for_similar_proteins,
    set_hmm_outdegree,
)
from rich.progress import Progress
from rich.progress import (
    Progress,
)
from math import ceil
import logging
import argparse

import pandas as pd
from socialgene.utils.logging import log
from textdistance import smith_waterman

logging.getLogger("neo4j").setLevel(logging.WARNING)
logging.getLogger().setLevel(logging.INFO)


parser = argparse.ArgumentParser(
    description="Create a clustermap.js file for an input BGC"
)

parser.add_argument(
    "--input",
    help="Output directory filepath",
)
parser.add_argument(
    "--output",
    help="Output directory filepath",
)
parser.add_argument(
    "--hmmdir",
    help="HMM file directory",
)


def main(input_tuple):
    try:
        main2(input_tuple)
        return 2
    except:
        with open(Path(input_tuple[1], f"{Path(input_tuple[0]).stem}.json"), "w") as fp:
            pass
        return 1


def main2(input_tuple):
    gbk_path = input_tuple[0]
    json_dir = input_tuple[1]
    hmm_dir = input_tuple[2]
    # hmm file with cutoffs
    h1 = Path(hmm_dir, "socialgene_nr_hmms_file_with_cutoffs_1_of_1.hmm")
    # hmm file without cutoffs
    h2 = Path(hmm_dir, "socialgene_nr_hmms_file_without_cutoffs_1_of_1.hmm")
    ########################################
    # Setup hmmsearch
    ########################################
    # hmmpress
    HMMER().hmmpress(h1)
    HMMER().hmmpress(h2)

    ########################################
    # Read input BGC
    ########################################

    sg_object = SocialGene()

    # parse input BGC
    sg_object.parse(gbk_path)

    # fragile implementation; add string input BGC id to keep track/separate input BGC from those added to class later
    input_bgc_name = list(sg_object.assemblies.keys())[0]
    modified_input_bgc_name = f"socialgene_query_{input_bgc_name}"

    # change input assembly name
    sg_object.assemblies[modified_input_bgc_name] = sg_object.assemblies.pop(
        input_bgc_name
    )
    sg_object.assemblies[modified_input_bgc_name].id = modified_input_bgc_name

    input_bgc_nuc_id = list(sg_object.assemblies[modified_input_bgc_name].loci.keys())[
        0
    ]

    # Annotate input BGC proteins with HMM models using external hmmscan
    sg_object.annotate(
        hmm_directory=hmm_dir,
        use_neo4j_precalc=True,
    )

    # Check if HMM Neo4j nodes have an outdegree property; if not, calculate and add to node
    if not check_for_hmm_outdegree():
        set_hmm_outdegree()

    ########################################

    # Rank input proteins by how many (:hmm)-[:ANNOTATES]->(:protein) relationships will have to be traversed
    # Don't reduce, but limit to 25 searched proteins per input BGC
    # input_protein_domain_df = prioritize_input_proteins(sg_object, reduce_by=1, max_out=50)

    protein_domain_dict = get_lowest_outdegree_model_per_protein(sg_object)

    searched_protein_set = set(protein_domain_dict.keys())
    sg_object.assemblies[modified_input_bgc_name].loci[
        input_bgc_nuc_id
    ].sort_by_middle()

    # input_sorted_protein_hashlist = [
    #     i.protein_hash
    #     for i in sg_object.assemblies[modified_input_bgc_name]
    #     .loci[input_bgc_nuc_id]
    #     .features
    # ]

    # input_proteinget_lowest_outdegree_model_per_protein

    input_protein_domain_df = protein_domain_dict

    # CF_002362315.1

    initial_search_results = search_for_similar_proteins(protein_domain_dict)

    #################################33
    # prioritized_nucleotide_df = count_matches(initial_search_results).sort_values(
    #     by=["locus_contains_n_queries", "assembly_contains_n_queries"],
    #     ascending=[False, False],
    #     na_position="first",
    # )

    # must_find_more_than = len(input_protein_domain_df) * 0.75

    # prioritized_nucleotide_df = prioritized_nucleotide_df[
    #     (prioritized_nucleotide_df["assembly_contains_n_queries"] >= must_find_more_than)
    # ]

    # if must_find_more_than > 3:
    #     nuc_gt = must_find_more_than * 0.75
    # else:
    #     nuc_gt = 2

    # prioritized_nucleotide_df = prioritized_nucleotide_df[
    #         (prioritized_nucleotide_df["locus_contains_n_queries"] >= nuc_gt)
    #     ]
    #################################33

    must_find_more_than_x_per_nuc_sequence = ceil(len(input_protein_domain_df) * 0.5)

    assemblies_must_have_x_matches = ceil(len(input_protein_domain_df)) * 0.75

    links_must_have_mod_score_above = 0.8

    #####################################
    # SECOND PASS
    #####################################

    ###################################################################################################################

    ###################################################################################################################
    assembly_counts = (
        initial_search_results.groupby("assembly_uid")["query_prot_uid"]
        .nunique()
        .reset_index()
    )

    reduced = initial_search_results.merge(
        assembly_counts[
            assembly_counts.query_prot_uid >= assemblies_must_have_x_matches
        ]["assembly_uid"],
        how="inner",
        on="assembly_uid",
    )

    reduced = reduced.sort_values(by=["n_start"], ascending=[True], na_position="first")
    reduced.reset_index(inplace=True, drop=True)

    # #####################################
    # # first pass
    # #####################################
    # big_clusters = []

    # with Progress(transient=True) as pg:
    #     task = pg.add_task("Finding clusters...", total=reduced["nucleotide_uid"].nunique())
    #     for name, group in reduced.groupby(["nucleotide_uid"]):
    #         reduced["cluster"] = (
    #             group["n_start"] > (group["n_end"].shift() + 10000)
    #         ).cumsum()
    #         reduced["unique_hits"] = group.groupby("cluster")["query_prot_uid"].transform(
    #             "nunique"
    #         )
    #         # big_clusters.append(group)
    #         pg.update(task, advance=1)

    ##################################################################################################################
    reduced["cluster"] = (
        reduced.groupby("nucleotide_uid")
        .apply(lambda group: group["n_start"] > (group["n_end"].shift() + 20000))
        .cumsum()
        .reset_index(level=[0])[0]
    )

    reduced["cluster_unique_hits"] = reduced.groupby(["nucleotide_uid", "cluster"])[
        "query_prot_uid"
    ].transform("nunique")

    reduced["assembly_unique_hits"] = reduced.groupby(["assembly_uid"])[
        "query_prot_uid"
    ].transform("nunique")

    reduced[reduced["assembly_uid"] == "GCF_002362315.1"]

    clusters_above_threshold = reduced[
        reduced["cluster_unique_hits"] > must_find_more_than_x_per_nuc_sequence
    ]
    clusters_above_threshold_final_df = (
        clusters_above_threshold.groupby(["assembly_uid", "nucleotide_uid", "cluster"])
        .agg({"n_start": "min", "n_end": "max"})
        .reset_index()
    )
    reduced.to_csv(Path(json_dir, f"{Path(gbk_path).stem}_reduced_1.tsv"), sep="\t")
    initial_search_results.to_csv(
        Path(json_dir, f"{Path(gbk_path).stem}_initial_search_results.tsv"), sep="\t"
    )
    clusters_above_threshold_final_df.to_csv(
        Path(json_dir, f"{Path(gbk_path).stem}_clusters_above_threshold_final_df.tsv"),
        sep="\t",
    )

    reduced.drop(clusters_above_threshold.index, inplace=True)
    # reduce to only genome assemblies that hada large clusters
    reduced = reduced[
        reduced["assembly_uid"].isin(clusters_above_threshold["assembly_uid"].unique())
    ]
    reduced = reduced[
        reduced["assembly_uid"].isin(clusters_above_threshold["assembly_uid"].unique())
    ]

    reduced["cluster_unique_hits"] = reduced.groupby(["nucleotide_uid", "cluster"])[
        "query_prot_uid"
    ].transform("nunique")

    reduced.to_csv(Path(json_dir, f"{Path(gbk_path).stem}_reduced_2.tsv"), sep="\t")

    # small_groups = []
    # for name, group in small_clusters.groupby(["assembly_uid"]):
    #     lg = group[group["new_hits"] == group["new_hits"].max()]
    #     lg = (
    #         lg.groupby("nucleotide_uid")
    #         .agg({"n_start": "min", "n_end": "max"})
    #         .reset_index()
    #     )
    #     small_groups.append(
    #         lg.iloc[(lg["n_end"] - lg["n_start"]).idxmin()]
    #         .to_frame()
    #         .transpose()
    #         .reset_index(drop=True)
    #     )
    #     # reduced.drop(lg.index, inplace=True)

    # small_groups = pd.concat(small_groups)

    #######################################
    # By here we have found all genomes with initial BGCs
    #######################################

    # big_clusters2 = (
    #     main_clusters.groupby("nucleotide_uid")
    #     .agg({"n_start": "min", "n_end": "max"})
    #     .reset_index()
    # )

    # sm_clusters2 = small_groups

    # Add big clusters

    with Progress(transient=True) as pg:
        task = pg.add_task(
            "Adding best hits...", total=len(clusters_above_threshold_final_df)
        )
        for index, result in clusters_above_threshold_final_df.filter(
            ["nucleotide_uid", "n_start", "n_end"]
        ).iterrows():
            sg_object.fill_given_locus_range(
                locus_uid=result[0], start=result[1] - 5000, end=result[2] + 5000
            )
            pg.update(task, advance=1)

    # with Progress(transient=True) as pg:
    #     task = pg.add_task("Filling proteins found elsewhere...", total=len(other_proteins))
    #     for index, result in other_proteins.filter(['nucleotide_uid','n_start','n_end']).iterrows():
    #         sg_object.fill_given_locus_range(
    #             locus_uid=result[0], start=result[1], end=result[2]
    #         )
    #         pg.update(task, advance=1)

    #######################################

    # Add small clusters

    # with Progress(transient=True) as pg:
    #     task = pg.add_task("Adding best hits...", total=len(sm_clusters2))
    #     for index, result in sm_clusters2.filter(
    #         ["nucleotide_uid", "n_start", "n_end"]
    #     ).iterrows():
    #         sg_object.fill_given_locus_range(
    #             locus_uid=result[0], start=result[1] - 5000, end=result[2] + 5000
    #         )
    #         pg.update(task, advance=1)

    sg_object.protein_comparison = []

    query_assembly = [
        v for k, v in sg_object.assemblies.items() if k.startswith("socialgene_query_")
    ][0]
    query_proteins = query_assembly.get_all_proteins()

    target_proteins = set()
    for k, v in sg_object.assemblies.items():
        if not k.startswith("socialgene_query_"):
            target_proteins.update(v.get_all_proteins())

    sg_object.bro(queries=query_proteins, targets=target_proteins, append=True)

    #######################################

    sg_object.protein_comparison_to_df()
    # sg_object.protein_comparison = sg_object.protein_comparison[
    #     sg_object.protein_comparison.mod_score > 1.4
    # ]

    sg_object.protein_comparison = sg_object.protein_comparison[
        sg_object.protein_comparison.mod_score > links_must_have_mod_score_above
    ]

    groupdict = (
        sg_object.protein_comparison.groupby("query")["target"].apply(list).to_dict()
    )

    group_dict_info = {
        f.protein_hash: (f.protein_id, f.description)
        for f in sg_object.assemblies[modified_input_bgc_name]
        .loci[input_bgc_nuc_id]
        .features
    }

    ######################## ORDER THE OUTPUT BGCs
    # get {nucleotide_hash: assembly_uid}
    zz = {
        sg_object._create_internal_locus_id(k, k2): k
        for k, v in sg_object.assemblies.items()
        for k2 in v.loci.keys()
    }

    clusters_above_threshold.sort_values("cluster_unique_hits", ascending=False)[
        "assembly_uid"
    ].unique()

    # a = clusters_above_threshold.groupby(["assembly_uid", "nucleotide_uid", "cluster"])[
    #     "query_prot_uid"
    # ].apply(lambda x: mod_score(list(x), input_sorted_protein_hashlist)["mod_score"])

    # a = a.reset_index()
    # a = a.sort_values("query_prot_uid", ascending=True)

    # a = a["assembly_uid"].unique()
    a = clusters_above_threshold.sort_values("cluster_unique_hits", ascending=False)[
        "assembly_uid"
    ].unique()

    order2 = [modified_input_bgc_name] + list(a)

    ########################

    ######################## Modify Assembly names

    for k, v in sg_object.assemblies.items():
        v.fill_taxonomy_from_db()
        v.fill_properties()

    for k, v in sg_object.assemblies.items():
        if v.taxonomy.genus_:
            v.name = f"{v.id}_{v.taxonomy.genus_}_{v.info['culture_collection']}"
        else:
            v.name = f"{v.id}_{v.info['culture_collection']}"
        try:
            a_hits = reduced[reduced["assembly_uid"] == k]["assembly_unique_hits"].iloc[
                0
            ]
            a_perc = round(a_hits / len(searched_protein_set) * 100, 0)
            n_hits = (
                clusters_above_threshold[clusters_above_threshold["assembly_uid"] == k]
                .filter(["cluster", "cluster_unique_hits"])
                .drop_duplicates()["cluster_unique_hits"]
                .sum()
            )
            n_perc = round(n_hits / len(searched_protein_set) * 100, 0)
            v.name = f"{v.name}_{str(n_perc)}%_{str(a_perc)}%"
        except:
            pass

    cmap = Clustermap()

    cmap.write(
        sg_object,
        groupdict=groupdict,
        group_dict_info=group_dict_info,
        assembly_order=order2,
        outpath=Path(json_dir, f"{Path(gbk_path).stem}.json"),
    )

    sg_object.ferment_pickle(Path(json_dir, f"{Path(gbk_path).stem}.pickle"))


if __name__ == "__main__":
    args = parser.parse_args()
    gb_dir = args.input
    json_dir = args.output
    hmm_dir = args.hmmdir

    gb_files = [
        i
        for i in Path(gb_dir).glob("*.gbk")
        if not Path(json_dir, f"{i.stem}.json").exists()
    ]

    p = Pool(40)

    for i in p.imap_unordered(main, itertools.product(gb_files, [json_dir], [hmm_dir])):
        pass
