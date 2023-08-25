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


mibig_id = "BGC0000001"
gbk_path = f"/home/chase/Documents/data/mibig/3_1/mibig_gbk_3.1/{mibig_id}.gbk"

hmm_dir = (
    "/home/chase/Documents/socialgene_data/streptomyces/socialgene_per_run/hmm_cache"
)

# hmm_dir='/home/chase/Downloads/ttt/hmm'


# hmm file with cutoffs
h1 = Path(hmm_dir, "socialgene_nr_hmms_file_with_cutoffs_1_of_1.hmm.gz")
# hmm file without cutoffs
h2 = Path(hmm_dir, "socialgene_nr_hmms_file_without_cutoffs_1_of_1.hmm.gz")

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
sg_object.assemblies[f"socialgene_query_{input_bgc_name}"] = sg_object.assemblies.pop(input_bgc_name)
sg_object.assemblies[f"socialgene_query_{input_bgc_name}"].id=f"socialgene_query_{input_bgc_name}"

# Annotate input BGC proteins with HMM models using external hmmscan
sg_object.annotate(
    hmm_directory=hmm_dir,
    use_neo4j_precalc=True,
)

########################################

########################################
if not check_for_hmm_outdegree():
    set_hmm_outdegree()


input_protein_domain_df = prioritize_input_proteins(sg_object, reduce_by=1, max_out=50)

initial_search_results = search_for_similar_proteins(input_protein_domain_df, sg_object)





prioritized_nucleotide_seqs = count_matches_per_nucleotide_sequence(
    initial_search_results
)


prioritized_nucleotide_seqs = prioritized_nucleotide_seqs[
    prioritized_nucleotide_seqs.count_of_matched_inputs
    >= len(input_protein_domain_df.query_prot_uid) * 0.5
]

bro = []

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
            bro.append(result)
        pg.update(task, advance=1)


bro = pd.concat(bro)

a = count_matches_per_nucleotide_sequence(bro)
a = a[a.count_of_matched_inputs > 1]

bro2 = []

with Progress(transient=True) as pg:
    task = pg.add_task("Processing initial matches...", total=len(a))
    for name, group in bro[bro["nuc_uid"].isin(list(a.nucleotide_uid))].groupby(
        ["nuc_uid"]
    ):
        result = prioritize_cluster_windows(group)
        if result:
            bro2.append(result)
        pg.update(task, advance=1)


groupdict = bro.groupby('query_prot_uid')['target_prot_uid'].apply(list).to_dict()

sg_res_object = sg_object


with Progress(transient=True) as pg:
    task = pg.add_task("Progress...", total=len(bro2))
    for result in bro2:
        sg_object.fill_given_locus_range(
            locus_uid=result[0], start=result[1] - 100, end=result[2] + 100
        )
        pg.update(task, advance=1)

sg_res_object.protein_comparison = []
sg_res_object.compare_proteins(append=True, cpus=20)
sg_res_object.protein_comparison_to_df()
sg_res_object.protein_comparison = sg_res_object.protein_comparison[
    sg_res_object.protein_comparison.mod_score >1.4
]


import json



a=Clustermap()

z=a.doit(sg_object, groupdict=groupdict, group_dict_info=group_dict_info, assembly_order=list(sg_object.assemblies.keys()))


with open("/home/chase/Downloads/ttt/tempppp/clinker/clinker/plot/data.json", "w") as outfile:
    json.dump(z, outfile)







########################################################################################################################
############################################################
############################################################
############################################################
cm_object = Clustermap()
cm_object.create_clustermap_uuids(sg_object=sg_res_object)
cm_object.add_cluster(sg_object=sg_res_object)
cm_object.add_groups(sg_object=sg_res_object, cutoff=0.8)
cm_object.add_links(sg_object=sg_res_object)
cm_object.write_clustermap("/home/chase/Downloads/ttt/data.json")


input_group_data = (
    initial_search_results.groupby("query_prot_uid")["target_prot_uid"]
    .apply(set)
    .reset_index()
)


input_group_data = {
    row["query_prot_uid"]: list(row["target_prot_uid"])
    for index, row in input_group_data.iterrows()
}

input_group_data

AXO35216.1
