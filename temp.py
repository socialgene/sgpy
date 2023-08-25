from socialgene.clustermap.clustermap import Clustermap
from socialgene.hmm.hmmer import HMMER
from socialgene.base.socialgene import SocialGene
from pathlib import Path
import pandas as pd

from socialgene.scoring.search import (
    check_for_hmm_outdegree,
    count_matches_per_nucleotide_sequence,
    prioritize_cluster_windows,
    prioritize_input_proteins,
    search_for_similar_proteins,
    set_hmm_outdegree,
)
from rich.progress import Progress

mibig_id = "BGC0001850"
gbk_path = f"/home/chase/Documents/data/mibig/3_1/mibig_gbk_3.1/{mibig_id}.gbk"

# hmm_dir = (
#     "/home/chase/Documents/socialgene_data/streptomyces/socialgene_per_run/hmm_cache"
# )

hmm_dir = "/home/chase/Downloads/ttt/hmm"


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

########################################

########################################
if not check_for_hmm_outdegree():
    set_hmm_outdegree()


input_protein_domain_df = prioritize_input_proteins(sg_object, reduce_by=1, max_out=100)

initial_search_results = search_for_similar_proteins(input_protein_domain_df, sg_object)


prioritized_nucleotide_seqs = count_matches_per_nucleotide_sequence(
    initial_search_results
)


prioritized_nucleotide_seqs = prioritized_nucleotide_seqs[
    prioritized_nucleotide_seqs.count_of_matched_inputs
    >= len(input_protein_domain_df) * 0.5
]

intermediate_df_1 = []

with Progress(transient=True) as pg:
    task = pg.add_task(
        "Processing initial matches...", total=len(prioritized_nucleotide_seqs)
    )
    for name, group in initial_search_results[
        initial_search_results["nuc_uid"].isin(
            list(prioritized_nucleotide_seqs.nucleotide_uid)
        )
    ].groupby(["nuc_uid"]):
        result = prioritize_cluster_windows(group, return_df_not_tuple=True)
        if isinstance(result, pd.DataFrame):
            intermediate_df_1.append(result)
        pg.update(task, advance=1)


intermediate_df_1 = pd.concat(intermediate_df_1)

df_nucuid_n_matches = count_matches_per_nucleotide_sequence(intermediate_df_1)
df_nucuid_n_matches = df_nucuid_n_matches[
    df_nucuid_n_matches.count_of_matched_inputs > len(input_protein_domain_df) * 0.5
]

intermediate_df_2 = []

with Progress(transient=True) as pg:
    task = pg.add_task("Processing initial matches...", total=len(df_nucuid_n_matches))
    for name, group in intermediate_df_1[
        intermediate_df_1["nuc_uid"].isin(list(df_nucuid_n_matches.nucleotide_uid))
    ].groupby(["nuc_uid"]):
        result = prioritize_cluster_windows(group)
        if result:
            intermediate_df_2.append(result)
        pg.update(task, advance=1)


groupdict = (
    intermediate_df_1.groupby("query_prot_uid")["target_prot_uid"].apply(list).to_dict()
)


with Progress(transient=True) as pg:
    task = pg.add_task("Progress...", total=len(intermediate_df_2))
    for result in intermediate_df_2:
        sg_object.fill_given_locus_range(
            locus_uid=result[0], start=result[1] - 10000, end=result[2] + 10000
        )
        pg.update(task, advance=1)

sg_object.protein_comparison = []
sg_object.compare_proteins(append=True, cpus=20)
sg_object.protein_comparison_to_df()
sg_object.protein_comparison = sg_object.protein_comparison[
    sg_object.protein_comparison.mod_score > 1.4
]


group_dict_info = {
    f.protein_hash: (f.protein_id, f.description)
    for f in sg_object.assemblies[modified_input_bgc_name]
    .loci[input_bgc_nuc_id]
    .features
}


# get {nucleotide_hash: assembly_uid}
zz = {
    sg_object._create_internal_locus_id(k, k2): k
    for k, v in sg_object.assemblies.items()
    for k2 in v.loci.keys()
}

order1 = [zz[i] for i in list(df_nucuid_n_matches.nucleotide_uid) if i in zz]

order2 = [modified_input_bgc_name] + order1


cmap = Clustermap()

cmap.write(
    sg_object,
    groupdict=groupdict,
    group_dict_info=group_dict_info,
    assembly_order=order2,
    outpath="/home/chase/Downloads/ttt/tempppp/clinker/clinker/plot/data.json",
)
