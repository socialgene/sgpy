from socialgene.clustermap.clustermap import Clustermap
from socialgene.hmm.hmmer import HMMER
from socialgene.base.socialgene import SocialGene
from pathlib import Path
import pandas as pd
from socialgene.neo4j.neo4j import GraphDriver

from socialgene.scoring.search import (
    check_for_hmm_outdegree,
    count_matches,
    count_matches_per_nucleotide_sequence,
    prioritize_cluster_windows,
    prioritize_input_proteins,
    get_lowest_outdegree_model_per_protein,
    search_for_similar_proteins,
    set_hmm_outdegree,
)
from rich.progress import Progress
from datetime import datetime
from socialgene.neo4j.neo4j import GraphDriver
from rich.progress import (
    Progress,
)
from math import ceil
import logging
import pandas as pd
from socialgene.utils.logging import log

logging.getLogger("neo4j").setLevel(logging.WARNING)
logging.getLogger().setLevel(logging.INFO)

mibig_id = "BGC0001646"
gbk_path = f"/home/chase/mibig_gbk/mibig_gbk_3.1/{mibig_id}.gbk"
hmm_dir = "/media/bigdrive2/chase/socialgene_v0.2.3/refseq/socialgene_per_run/hmm_cache"


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
sg_object.assemblies[modified_input_bgc_name] = sg_object.assemblies.pop(input_bgc_name)
sg_object.assemblies[modified_input_bgc_name].id = modified_input_bgc_name

input_bgc_nuc_id = list(sg_object.assemblies[modified_input_bgc_name].loci.keys())[0]

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


zz = assembly_counts[assembly_counts.query_prot_uid > assemblies_must_have_x_matches]

reduced = initial_search_results.merge(
    assembly_counts[assembly_counts.query_prot_uid > assemblies_must_have_x_matches][
        "assembly_uid"
    ],
    how="inner",
    on="assembly_uid",
)


reduced = reduced.sort_values(by=["n_start"], ascending=[True], na_position="first")
reduced.reset_index(inplace=True, drop=True)


#####################################
# first pass
#####################################
big_clusters = []
with Progress(transient=True) as pg:
    task = pg.add_task("Finding clusters...", total=reduced["nucleotide_uid"].nunique())
    for name, group in reduced.groupby(["nucleotide_uid"]):
        group["cluster"] = (
            group["n_start"] > (group["n_end"].shift() + 10000)
        ).cumsum()
        group["unique_hits"] = group.groupby("cluster")["query_prot_uid"].transform(
            "nunique"
        )
        big_clusters.append(group)
        pg.update(task, advance=1)


big_clusters = pd.concat(big_clusters)

main_clusters = big_clusters[
    big_clusters["unique_hits"] > must_find_more_than_x_per_nuc_sequence
]

reduced.drop(main_clusters.index, inplace=True)

reduced = reduced[
    reduced["assembly_uid"].isin(main_clusters["assembly_uid"])
].reset_index(drop=True)

reduced = reduced.sort_values(by=["n_start"], ascending=[True], na_position="first")
reduced.reset_index(inplace=True, drop=True)
reduced

#####################################
# second pass
#####################################
small_clusters = []
with Progress(transient=True) as pg:
    task = pg.add_task("Finding clusters...", total=reduced["assembly_uid"].nunique())
    for name, group in reduced.groupby(["assembly_uid"]):
        already_found_proteins = main_clusters[
            main_clusters["assembly_uid"] == name[0]
        ]["query_prot_uid"].unique()
        group["cluster"] = (
            group["n_start"] > (group["n_end"].shift() + 100000)
        ).cumsum()
        group["new_hits"] = group.groupby("cluster")["query_prot_uid"].transform(
            lambda x: len(set(x) - set(already_found_proteins))
        )
        small_clusters.append(group)
        pg.update(task, advance=1)


small_clusters = pd.concat(small_clusters)


small_groups = []
for name, group in small_clusters.groupby(["assembly_uid"]):
    lg = group[group["new_hits"] == group["new_hits"].max()]
    lg = (
        lg.groupby("nucleotide_uid")
        .agg({"n_start": "min", "n_end": "max"})
        .reset_index()
    )
    small_groups.append(
        lg.iloc[(lg["n_end"] - lg["n_start"]).idxmin()]
        .to_frame()
        .transpose()
        .reset_index(drop=True)
    )

    small_groups.append(lg)
    # reduced.drop(lg.index, inplace=True)


small_groups = pd.concat(small_groups)


#######################################
# By here we have found all genomes with initial BGCs
#######################################

big_clusters2 = (
    main_groups.groupby("nucleotide_uid")
    .agg({"n_start": "min", "n_end": "max"})
    .reset_index()
)

sm_clusters2 = small_groups


big_clusters2 = big_clusters2[0:20]

# Add big clusters

with Progress(transient=True) as pg:
    task = pg.add_task("Adding best hits...", total=len(big_clusters2))
    for index, result in big_clusters2.filter(
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
sg_object.protein_comparison_to_df()
# sg_object.protein_comparison = sg_object.protein_comparison[
#     sg_object.protein_comparison.mod_score > 1.4
# ]

sg_object.protein_comparison = sg_object.protein_comparison[
    sg_object.protein_comparison.mod_score > 0.8
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


df_nucuid_n_matches = (
    big_clusters.groupby(["nucleotide_uid"])
    .first()
    .filter(["nucleotide_uid", "unique_queries"])
    .reset_index(drop=False)
)

order1 = [zz[i] for i in list(df_nucuid_n_matches.nucleotide_uid) if i in zz]

order2 = [modified_input_bgc_name] + list(set(order1))
########################

######################## Modify Assembly names


for k, v in sg_object.assemblies.items():
    v.fill_taxonomy_from_db()
    v.fill_properties()

for k, v in sg_object.assemblies.items():
    if v.taxonomy.genus_:
        v.name = f"{v.id}__{v.taxonomy.genus_}__{v.info['culture_collection']}"
    else:
        v.name = f"{v.id}____{v.info['culture_collection']}"


cmap = Clustermap()

cmap.write(
    sg_object,
    groupdict=groupdict,
    group_dict_info=group_dict_info,
    assembly_order=order2,
    outpath="/home/chase/data.json",
)
