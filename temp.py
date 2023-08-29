import itertools
from socialgene.clustermap.clustermap import Clustermap
from socialgene.hmm.hmmer import HMMER
from socialgene.base.socialgene import SocialGene
from pathlib import Path
import pandas as pd
import argparse
from multiprocessing import Pool
from pathlib import Path
from socialgene.scoring.search import (
    check_for_hmm_outdegree,
    count_matches_per_nucleotide_sequence,
    prioritize_cluster_windows,
    prioritize_input_proteins,
    search_for_similar_proteins,
    set_hmm_outdegree,
)
from rich.progress import Progress


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

# gbk_path = f"/home/chase/mibig_gbk/mibig_gbk_3.1/{mibig_id}.gbk"
# hmm_dir = "/media/bigdrive2/chase/socialgene_v0.2.3/refseq/socialgene_per_run/hmm_cache"


def main(a):
    try:
        main2(a)
    except:
        pass
              

def main2(input_tuple):
    print(input_tuple)
    gbk_path = input_tuple[0]
    json_dir= input_tuple[1]
    hmm_dir = input_tuple[2] 
    
    json_path = Path(json_dir, f"{Path(gbk_path).stem}.json")
    if json_path.exists():
        return
    
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


    input_protein_domain_df = prioritize_input_proteins(sg_object, reduce_by=1.5, max_out=25)

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
            initial_search_results["nucleotide_uid"].isin(
                list(prioritized_nucleotide_seqs.nucleotide_uid)
            )
        ].groupby(["nucleotide_uid"]):
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
            intermediate_df_1["nucleotide_uid"].isin(list(df_nucuid_n_matches.nucleotide_uid))
        ].groupby(["nucleotide_uid"]):
            result = prioritize_cluster_windows(group)
            if result:
                intermediate_df_2.append(result)
            pg.update(task, advance=1)



    with Progress(transient=True) as pg:
        task = pg.add_task("Progress...", total=len(intermediate_df_2))
        for result in intermediate_df_2:
            sg_object.fill_given_locus_range(
                locus_uid=result[0], start=result[1] - 10000, end=result[2] + 10000
            )
            pg.update(task, advance=1)

    sg_object.protein_comparison = []


    query_assembly=[v for k,v in sg_object.assemblies.items() if k.startswith("socialgene_query_")][0]
    query_proteins = query_assembly.get_all_proteins()

    target_proteins=set()

    for k,v in sg_object.assemblies.items():
        if not k.startswith("socialgene_query_"):
           target_proteins.update(v.get_all_proteins())


    sg_object.bro(queries=query_proteins, targets=target_proteins, append=True)
    sg_object.protein_comparison_to_df()

    sg_object.protein_comparison = sg_object.protein_comparison[
        sg_object.protein_comparison.jaccard > 0.7
    ]


    groupdict = sg_object.protein_comparison.groupby("query")['target'].apply(list).to_dict()


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

    order1 = [zz[i] for i in list(df_nucuid_n_matches.nucleotide_uid) if i in zz]

    order2 = [modified_input_bgc_name] + order1
    ########################

    ######################## Modify Assembly names



    for k,v in sg_object.assemblies.items():
        v.fill_taxonomy_from_db()
        v.fill_properties()

    for k,v in sg_object.assemblies.items():
        if v.taxonomy.genus_:
            v.name=f"{v.id}__{v.taxonomy.genus_}__{v.info['culture_collection']}"
        else:
            v.name=f"{v.id}____{v.info['culture_collection']}"

    cmap = Clustermap()

    cmap.write(
        sg_object,
        groupdict=groupdict,
        group_dict_info=group_dict_info,
        assembly_order=order2,
        outpath=json_path,
    )

if __name__ == "__main__":
    args = parser.parse_args()
    gb_dir = args.input
    json_dir=args.output
    hmm_dir = args.hmmdir 
    p = Pool(60)
    p.map(main, itertools.product(Path(gb_dir).glob("*.gbk"), [json_dir], [hmm_dir]))






    