import multiprocessing
import pickle
from multiprocessing import Pool
from pathlib import Path

from socialgene.cli.search.sea import search_bgc
from socialgene.config import env_vars

env_vars["NEO4J_URI"] = "bolt://localhost:7687"


def run_individual_search(gbk_path):
    outpath_clinker = Path(
        "/home/chase/Downloads/clinker/clinker/plot/jsons", gbk_path.stem + ".json"
    )
    if not outpath_clinker.exists():
        try:
            a = search_bgc(
                input=gbk_path,
                hmm_dir="/home/chase/Documents/socialgene_data/v0_4_1/refseq_antismash_bgcs/socialgene_per_run/hmm_cache",
                outpath_clinker=outpath_clinker,
                use_neo4j_precalc=True,
                assemblies_must_have_x_matches=0.4,
                nucleotide_sequences_must_have_x_matches=0.4,
                gene_clusters_must_have_x_matches=0.4,
                break_bgc_on_gap_of=10000,
                target_bgc_padding=20000,
                max_domains_per_protein=3,
                max_outdegree=1000000,
                max_query_proteins=5,
                scatter=False,
                locus_tag_bypass_list=None,
                protein_id_bypass_list=None,
                only_culture_collection=True,
                frac=0.75,
                run_async=True,
                analyze_with="blastp",
            )
            with open(
                Path(
                    "/home/chase/Downloads/clinker/clinker/plot/pickles",
                    gbk_path.stem + ".pickle",
                ),
                "wb",
            ) as handle:
                pickle.dump(a, handle, protocol=pickle.HIGHEST_PROTOCOL)
        except Exception as e:
            with open(
                Path(
                    "/home/chase/Downloads/clinker/clinker/plot/errors",
                    gbk_path.stem + ".error",
                ),
                "w",
            ) as handle:
                handle.write(str(e))


with Pool(processes=multiprocessing.cpu_count() - 2) as pool:
    pool.map(
        run_individual_search,
        Path("/home/chase/Documents/data/mibig/3_1/mibig_gbk_3.1").rglob("*.gbk"),
    )
