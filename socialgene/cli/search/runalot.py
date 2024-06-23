#!/usr/bin/env python3

import pickle
from pathlib import Path
from socialgene.cli.search.sea import search_bgc
from socialgene.config import env_vars
import sys

env_vars["NEO4J_URI"] = "bolt://localhost:7687"

outdir = "/media/socialgene_nvme/culture_search_results"

Path(outdir).mkdir(exist_ok=True)
Path(outdir, "jsons").mkdir(exist_ok=True)
Path(outdir, "pickles").mkdir(exist_ok=True)
Path(outdir, "errors").mkdir(exist_ok=True)


def run_individual_search(gbk_path):
    outpath_clinker = Path(
        "/media/socialgene_nvme/culture_search_results/jsons", gbk_path.stem + ".json"
    )
    if not outpath_clinker.exists() and not Path(
                    "/media/socialgene_nvme/culture_search_results/errors",
                    gbk_path.stem + ".error",
                ).exists():
        try:
            a = search_bgc(
                input=gbk_path,
                hmm_dir="/media/socialgene_nvme/v0.4.1/refseq/socialgene_per_run/hmm_cache",
                outpath_clinker=outpath_clinker,
                use_neo4j_precalc=True,
                assemblies_must_have_x_matches=0.4,
                nucleotide_sequences_must_have_x_matches=0.4,
                gene_clusters_must_have_x_matches=0.4,
                break_bgc_on_gap_of=10000,
                target_bgc_padding=10000,
                max_domains_per_protein=3,
                max_outdegree=300000,
                max_query_proteins=10,
                scatter=True,
                locus_tag_bypass_list=None,
                protein_id_bypass_list=None,
                only_culture_collection=True,
                frac=0.75,
                run_async=True,
                analyze_with="blastp",
                blast_speed='ultra-sensitive',
            )
            with open(
                Path(
                    "/media/socialgene_nvme/culture_search_results/pickles",
                    gbk_path.stem + ".pickle",
                ),
                "wb",
            ) as handle:
                pickle.dump(a, handle, protocol=pickle.HIGHEST_PROTOCOL)
        except Exception as e:
            with open(
                Path(
                    "/media/socialgene_nvme/culture_search_results/errors",
                    gbk_path.stem + ".error",
                ),
                "w",
            ) as handle:
                handle.write(str(e))

if __name__ == "__main__":
    # run run_individual_search() using a command line input as the gbk_path
    gbk_path = Path(sys.argv[1])
    run_individual_search(gbk_path)
