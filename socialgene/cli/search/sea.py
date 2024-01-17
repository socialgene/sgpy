from typing import List

from socialgene.clustermap.serialize import SerializeToClustermap
from socialgene.config import env_vars
from socialgene.search.hmmer import SearchDomains
from socialgene.utils.logging import log

env_vars["NEO4J_URI"] = "bolt://localhost:7687"


def search_bgc(
    input,
    hmm_dir: str = None,
    use_neo4j_precalc: bool = True,
    assemblies_must_have_x_matches: float = 0.9,
    nucleotide_sequences_must_have_x_matches: float = 0.9,
    gene_clusters_must_have_x_matches: float = 0.9,
    break_bgc_on_gap_of: int = 10000,
    target_bgc_padding: int = 20000,
    max_domains_per_protein: float = None,
    max_outdegree: int = None,
    max_query_proteins: int = 5,
    scatter: bool = False,
    locus_tag_bypass_list: List[str] = None,
    protein_id_bypass_list: List[str] = None,
    only_culture_collection: bool = False,
    frac: float = 0.75,
    run_async: bool = False,
    analyze_with: str = "hmmer",
    outpath_clinker: str = None,
):
    log.info(f"Running search with args: {locals()}")
    search_object = SearchDomains(
        input=input,
        hmm_dir=hmm_dir,
        use_neo4j_precalc=use_neo4j_precalc,
        assemblies_must_have_x_matches=assemblies_must_have_x_matches,
        nucleotide_sequences_must_have_x_matches=nucleotide_sequences_must_have_x_matches,
        gene_clusters_must_have_x_matches=gene_clusters_must_have_x_matches,
        break_bgc_on_gap_of=break_bgc_on_gap_of,
        target_bgc_padding=target_bgc_padding,
    )
    search_object.outdegree_table
    search_object.prioritize_input_proteins(
        max_domains_per_protein=max_domains_per_protein,
        max_outdegree=max_outdegree,
        max_query_proteins=max_query_proteins,
        scatter=scatter,
        locus_tag_bypass_list=locus_tag_bypass_list,
        protein_id_bypass_list=protein_id_bypass_list,
    )
    search_object.outdegree_table
    search_object.search(
        only_culture_collection=only_culture_collection, frac=frac, run_async=run_async
    )
    # filters assemblies and nucleotide seqs that fall below threshold % of query proteins
    search_object.filter()
    # labels clusters with a unique id based on break_bgc_on_gap_of
    search_object.label_clusters()
    df = search_object._primary_bgc_regions()
    df["n_start"] = df["n_start"] - search_object.target_bgc_padding
    df["n_end"] = df["n_end"] + search_object.target_bgc_padding
    search_object._bgc_regions_to_sg_object(df)
    _ = search_object.sg_object.annotate_proteins_with_neo4j(
        protein_uids=None, annotate_all=True, progress=False
    )
    ########
    # add bgcs as gene_clusters to the locus objects
    for i, row in df.iterrows():
        search_object.sg_object.assemblies[row["assembly_uid"]].get_locus_by_uid(
            row["nucleotide_uid"]
        ).add_bgcs_by_start_end(
            start=row["n_start"] - search_object.target_bgc_padding,
            end=row["n_end"] + search_object.target_bgc_padding,
            uid=row["cluster"],
        )
    # add input bgc as gene_cluster to the locus objects
    if analyze_with == "hmmer":
        search_object.sg_object.annotate_proteins_with_neo4j(annotate_all=True)
    else:
        search_object.sg_object.add_sequences_from_neo4j()
    search_object._create_links(
        tool=analyze_with, argstring="--fast --max-hsps 1", cpus=10
    )
    search_object._choose_group()
    return search_object
    search_object._rank_order_bgcs()
    assemblies = list(
        dict.fromkeys([i.parent.parent for i in search_object.sorted_bgcs])
    )
    assemblies = [search_object.input_assembly] + assemblies
    z = SerializeToClustermap(
        sg_object=search_object.sg_object,
        sorted_bgcs=[search_object.input_bgc] + search_object.sorted_bgcs,
        link_df=search_object.link_df,
        group_df=search_object.group_df,
    )
    z.write(outpath_clinker)
