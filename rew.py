from math import ceil
from typing import List
from socialgene.base.socialgene import SocialGene
from socialgene.hmm.hmmer import HMMER
from socialgene.scoring.search import (
    _find_similar_protein_multiple,
    check_for_hmm_outdegree,
    get_outdegree_per_hmm_per_protein,
    set_hmm_outdegree,
    _find_similar_protein_multiple,
)
from pathlib import Path
import pandas as pd
import asyncio

input_gbk_path = "/home/chase/Documents/data/mibig/3_1/mibig_gbk_3.1/BGC0001850.gbk"
hmm_dir = (
    "/home/chase/Documents/socialgene_data/streptomyces/socialgene_per_run/hmm_cache"
)
# hmm file with cutoffs
h1 = Path(hmm_dir, "socialgene_nr_hmms_file_with_cutoffs_1_of_1.hmm")
# hmm file without cutoffs
h2 = Path(hmm_dir, "socialgene_nr_hmms_file_without_cutoffs_1_of_1.hmm")
hmm_file_list = [h1, h2]


def prep_hmms(hmm_file_list: List):
    for i in hmm_file_list:
        HMMER().hmmpress(i)


def read_input_bgc(gbk_path):
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
    return sg_object


prep_hmms(hmm_file_list)
sg_object = read_input_bgc(input_gbk_path)
# Check if HMM Neo4j nodes have an outdegree property; if not, calculate and add to node
if not check_for_hmm_outdegree():
    set_hmm_outdegree()


# Annotate input BGC proteins with HMM models using external hmmscan
sg_object.annotate(
    hmm_directory=hmm_dir,
    use_neo4j_precalc=True,
)


########################################
# Rank input proteins by how many (:hmm)-[:ANNOTATES]->(:protein) relationships will have to be traversed
# Don't reduce, but limit to 25 searched proteins per input BGC
# input_protein_domain_df = prioritize_input_proteins(sg_object, reduce_by=1, max_out=50)
domain_outdegree_df = get_outdegree_per_hmm_per_protein(sg_object)
protein_domain_df = domain_outdegree_df.groupby("protein_uid")["hmm_uid"].apply(list)

initial_search_results = await _find_similar_protein_multiple(
    protein_domain_df.to_dict()
)
# doing it for other columns will increase the mem usage
initial_search_results["query_prot_uid"] = initial_search_results[
    "query_prot_uid"
].astype("category")


########################################
nucleotide_sequences_must_have_x_matches = ceil(len(protein_domain_df) * 0.5)
gene_clusters_must_have_x_matches = ceil(len(protein_domain_df) * 0.5)
assemblies_must_have_x_matches = ceil(len(protein_domain_df)) * 0.75
links_must_have_mod_score_above = 0.4
########################################
assembly_counts = (
    initial_search_results.groupby("assembly_uid")["query_prot_uid"]
    .nunique()
    .reset_index()
)
# filter, keeping assemblies with >= assemblies_must_have_x_matches
reduced = initial_search_results.merge(
    assembly_counts[assembly_counts.query_prot_uid >= assemblies_must_have_x_matches][
        "assembly_uid"
    ],
    how="inner",
    on="assembly_uid",
)

# filter, keeping nucleotide seqs with >= assemblies_must_have_x_matches
nuc_counts = reduced.groupby("nucleotide_uid")["query_prot_uid"].nunique().reset_index()
reduced = reduced.merge(
    nuc_counts[nuc_counts.query_prot_uid >= nucleotide_sequences_must_have_x_matches][
        "nucleotide_uid"
    ],
    how="inner",
    on="nucleotide_uid",
)
# sort by start position of genes
reduced = reduced.sort_values(by=["n_start"], ascending=[True], na_position="first")
reduced.reset_index(inplace=True, drop=True)
########################################
########################################
########################################

##################################################################################################################

########################################
# Generate gene clusters
########################################

MAX_GAP = 20000

# label clusters
reduced["cluster"] = (
    reduced.groupby("nucleotide_uid")
    .apply(lambda group: group["n_start"] > (group["n_end"].shift() + MAX_GAP))
    .cumsum()
    .reset_index(level=[0])[0]
)

# count the number of unique input/query protein hits within each cluster
reduced["cluster_unique_hits"] = reduced.groupby(["nucleotide_uid", "cluster"])[
    "query_prot_uid"
].transform("nunique")

# count the number of unique input/query protein hits within each genome/assembly
reduced["assembly_unique_hits"] = reduced.groupby(["assembly_uid"])[
    "query_prot_uid"
].transform("nunique")


########################################
# Prune gene clusters
########################################

clusters_above_threshold = reduced[
    reduced["cluster_unique_hits"] > gene_clusters_must_have_x_matches
]


clusters_above_threshold_final_df = (
    clusters_above_threshold.groupby(["assembly_uid", "nucleotide_uid", "cluster"])
    .agg({"n_start": "min", "n_end": "max"})
    .reset_index()
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

########################################
# Retreive gene cluster data from database
########################################

progress_bar = Progress(
    TextColumn(f"Pulling target cluster data from database..."),
    BarColumn(),
    MofNCompleteColumn(),
    TextColumn("• Time elapsed "),
    TimeElapsedColumn(),
    transient=True,
)
with progress_bar as pg:
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

# Drop likely cross-origin proteins
for ak, av in sg_object.assemblies.items():
    for nk, nv in av.loci.items():
        sg_object.assemblies[ak].loci[nk].features = set(
            [i for i in nv.features if abs(i.end - i.start) < 100000]
        )


########################################
# Compare proteins input_bgc proteins against those found on initial pass
########################################
# First compare input proteins to all others to ensure groups
a = CompareDomains()
query_proteins = sg_object.assemblies[modified_input_bgc_name].get_all_proteins()
a.compare_many_to_many(
    (v for k, v in sg_object.proteins.items() if k in query_proteins),
    (v for k, v in sg_object.proteins.items() if k not in query_proteins),
)
df = a.df
df = df[df.score > links_must_have_mod_score_above]

########################################
# Compare proteins all-vs-all
########################################
a = CompareDomains()
query_proteins = sg_object.assemblies[modified_input_bgc_name].get_all_proteins()
a.compare_all_to_all_parallel(
    (v for k, v in sg_object.proteins.items() if k not in query_proteins), 22
)
df2 = a.df
df2 = df2[df2.score > links_must_have_mod_score_above]

########################################
# Compare proteins all-vs-all
########################################

df.groupby("query")["target"].apply(list).reset_index(name="new")

df = df[df.score > links_must_have_mod_score_above]

groupdict = df.groupby("target")["query"].apply(list).to_dict()

group_dict_info = {
    f.protein_hash: (f.protein_id, f.description)
    for f in sg_object.assemblies[modified_input_bgc_name]
    .loci[input_bgc_nuc_id]
    .features
}

########################################
# ORDER THE OUTPUT BGCs
########################################
# get {nucleotide_hash: assembly_uid}
zz = {
    sg_object._create_internal_locus_id(k, k2): k
    for k, v in sg_object.assemblies.items()
    for k2 in v.loci.keys()
}

clusters_above_threshold.sort_values("cluster_unique_hits", ascending=False)[
    "assembly_uid"
].unique()

a = clusters_above_threshold.sort_values("cluster_unique_hits", ascending=False)[
    "assembly_uid"
].unique()

order2 = [modified_input_bgc_name] + list(a)
########################################
# Create and write clustermap json
########################################

for k, v in sg_object.assemblies.items():
    v.fill_taxonomy_from_db()
    v.fill_properties()

# Modify names of assemblies for their appearance in clustermap
for k, v in sg_object.assemblies.items():
    if v.taxonomy.genus_:
        v.name = f"{v.id}_{v.taxonomy.genus_}_{v.info['culture_collection']}"
    else:
        v.name = f"{v.id}_{v.info['culture_collection']}"
    try:
        a_hits = reduced[reduced["assembly_uid"] == k]["assembly_unique_hits"].iloc[0]
        a_perc = round(a_hits / len(searched_protein_hash_set) * 100, 0)
        n_hits = (
            clusters_above_threshold[clusters_above_threshold["assembly_uid"] == k]
            .filter(["cluster", "cluster_unique_hits"])
            .drop_duplicates()["cluster_unique_hits"]
            .sum()
        )
        n_perc = round(n_hits / len(searched_protein_hash_set) * 100, 0)
        v.name = f"{v.name}_{str(n_perc)}%_{str(a_perc)}%"
    except:
        pass

cmap = Clustermap()

cmap.write(
    sg_object,
    groupdict=groupdict,
    group_dict_info=group_dict_info,
    assembly_order=order2,
    outpath=Path("/home/chase/data.json"),
)
