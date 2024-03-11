import pandas as pd

from socialgene.base.socialgene import SocialGene
from socialgene.cli.search.sea import search_bgc
from socialgene.config import env_vars
from socialgene.neo4j.neo4j import GraphDriver
from socialgene.utils.logging import log

env_vars["NEO4J_URI"] = "bolt://localhost:7687"


with GraphDriver() as db:
    bgcs = db.run(
        """
        MATCH (n:nucleotide)
        WHERE n.external_id STARTS WITH "BGC"
        RETURN n.uid as uid
        """,
    ).value()


counter = 0
threshold = 0.7

for bgc in bgcs:
    counter += 1
    log.warning(f"Running {counter} of {len(bgcs)}")
    sg = SocialGene()
    sg.fill_given_locus_range(locus_uid=bgc, start=float("-inf"), end=float("inf"))
    a = search_bgc(
        input=sg,
        hmm_dir="/home/chase/Documents/socialgene_data/v0_4_1/refseq_antismash_bgcs/socialgene_per_run/hmm_cache",
        outpath_clinker=None,
        use_neo4j_precalc=True,
        assemblies_must_have_x_matches=threshold,
        nucleotide_sequences_must_have_x_matches=threshold,
        gene_clusters_must_have_x_matches=threshold,
        break_bgc_on_gap_of=20000,
        target_bgc_padding=20000,
        max_domains_per_protein=3,
        max_outdegree=500000,
        max_query_proteins=10,
        scatter=False,
        locus_tag_bypass_list=None,
        protein_id_bypass_list=None,
        only_culture_collection=False,
        frac=0.75,
        run_async=True,
        analyze_with="blastp",
    )
    temp = pd.merge(
        a._compare_bgcs_by_jaccard_and_levenshtein(),
        a._compare_bgcs_by_median_bitscore(),
        left_on="query_gene_cluster",
        right_on="target_gene_cluster",
        how="inner",
    ).sort_values(by=["modscore", "score"], ascending=False)

    temp["nuc"] = temp["query_gene_cluster_x"].apply(lambda x: x.parent.uid)
    # get best hit per nuc sequence
    temp = temp.groupby("nuc").head(1).reset_index(drop=True)
    temp = temp[["query_gene_cluster_x", "jaccard", "levenshtein", "modscore", "score"]]
    temp = temp[temp["jaccard"] > threshold]
    temp["start"] = temp["query_gene_cluster_x"].apply(
        lambda x: min(i.start for i in x.features)
    )
    temp["end"] = temp["query_gene_cluster_x"].apply(
        lambda x: max(i.end for i in x.features)
    )
    temp["query_gene_cluster_x"] = temp["query_gene_cluster_x"].apply(
        lambda x: x.parent.uid
    )
    # remove rows if query_gene_cluster_x is the same as bgc
    temp = temp[temp["query_gene_cluster_x"] != bgc]
    temp = temp.to_dict(orient="records")
    if temp:
        log.info(f"Linking {len(temp)} assemblies to {bgc}")
        with GraphDriver() as db:
            db.run(
                """
                MATCH (n:nucleotide)
                WHERE n.uid=$bgc_id
                UNWIND $to_link as to_link
                MATCH (m:nucleotide)
                WHERE m.uid=to_link['query_gene_cluster_x']
                MERGE (n)-[r:LINKS_TO]->(m)
                    ON CREATE SET r.jaccard=to_link['jaccard'], r.levenshtein=to_link['levenshtein'], r.modscore=to_link['modscore'], r.score=to_link['score'], r.start=to_link['start'], r.end=to_link['end']
                """,
                bgc_id=bgc,
                to_link=temp,
            )
