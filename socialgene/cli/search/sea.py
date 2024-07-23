from typing import List

from socialgene.clustermap.serialize import SerializeToClustermap
from socialgene.search.hmmer import SearchDomains
from socialgene.utils.logging import log


def search_bgc(
    input,
    hmm_dir: str = None,
    use_neo4j_precalc: bool = True,
    assemblies_must_have_x_matches: float = 0.6,
    nucleotide_sequences_must_have_x_matches: float = 0.6,
    gene_clusters_must_have_x_matches: float = 0.6,
    break_bgc_on_gap_of: int = 10000,
    target_bgc_padding: int = 20000,
    max_domains_per_protein: float = 3,
    max_outdegree: int = 100000,
    max_query_proteins: int = 10,
    scatter: bool = False,
    locus_tag_bypass_list: List[str] = None,
    protein_id_bypass_list: List[str] = None,
    only_culture_collection: bool = False,
    frac: float = 0.75,
    run_async: bool = True,
    analyze_with: str = "blastp",
    outpath_clinker: str = None,
    blast_speed: str = "fast",
    blast_threads: int = 1,
    adaptive_padding: bool = True,
):
    log.info(f"Running search with args: {locals()}")

    match blast_speed:
        case "fast":
            blastarg = "--fast"
        case "sensitive":
            blastarg = "--sensitive"
        case "ultra-sensitive":
            blastarg = "--ultra-sensitive"
        case _:
            raise ValueError(f"Invalid blast speed: {blast_speed}")

    search_object = SearchDomains(
        input=input,
        hmm_dir=hmm_dir,
        use_neo4j_precalc=use_neo4j_precalc,
        assemblies_must_have_x_matches=assemblies_must_have_x_matches,
        nucleotide_sequences_must_have_x_matches=nucleotide_sequences_must_have_x_matches,
        gene_clusters_must_have_x_matches=gene_clusters_must_have_x_matches,
        break_bgc_on_gap_of=break_bgc_on_gap_of,
        target_bgc_padding=target_bgc_padding,
        adaptive_padding=adaptive_padding,
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
    if search_object.raw_search_results_df.empty:
        log.warning("Stopping, no search results to process")
        raise
    search_object.filter()
    # labels clusters with a unique id based on break_bgc_on_gap_of
    search_object.label_clusters()
    df = search_object._primary_bgc_regions()
    log.info(
        f"First pass resulted in {df.assembly_uid.nunique()} assemblies, {df.nucleotide_uid.nunique()} nucleotide sequences had {df.cluster.nunique()} putative BGCs"
    )

    df["n_start"] = df["n_start"] - search_object.search_padding
    df["n_end"] = df["n_end"] + search_object.search_padding
    search_object._bgc_regions_to_sg_object(df)

    ########
    # add bgcs as gene_clusters to the locus objects
    for i, row in df.iterrows():
        search_object.sg_object.assemblies[row["assembly_uid"]].get_locus_by_uid(
            row["nucleotide_uid"]
        ).add_bgcs_by_start_end(
            start=row["n_start"] - search_object.search_padding,
            end=row["n_end"] + search_object.search_padding,
            uid=row["cluster"],
        )
    # add input bgc as gene_cluster to the locus objects
    if analyze_with == "hmmer":
        _ = search_object.sg_object.annotate_proteins_with_neo4j(
            protein_uids=None, annotate_all=True, progress=False
        )
    else:
        search_object.sg_object.add_sequences_from_neo4j()

    # Find RBH between input BGC and putative BGCs
    if analyze_with == "blastp":
        argstring = f"{blastarg} --max-hsps 1"
    else:
        argstring = ""
    search_object._create_links(
        tool=analyze_with, argstring=argstring, threads=blast_threads
    )
    # Assign protein groups for the clinker plot legend
    search_object._choose_group()

    if not outpath_clinker:
        return search_object
    # return search_object
    if gene_clusters_must_have_x_matches > 1:
        # scale to 0-1
        gene_clusters_must_have_x_matches = gene_clusters_must_have_x_matches / len(
            search_object.input_bgc.proteins
        )
    for group, df in search_object.link_df.groupby("target_gene_cluster"):
        minstart = (
            df.target_feature.apply(lambda x: x.start).min()
            - search_object.target_bgc_padding
        )
        group.features = {i for i in group.features if i.start > minstart}
        maxend = (
            df.target_feature.apply(lambda x: x.end).max()
            + search_object.target_bgc_padding
        )
        group.features = {i for i in group.features if i.end < maxend}

    assemblies = search_object._rank_order_bgcs(
        threshold=gene_clusters_must_have_x_matches
    )
    assemblies = [i.parent.parent for i in assemblies]
    ids = [i.uid for i in assemblies]
    # log.warning(ids)
    if len(ids) == 1 and ids[0] == search_object.input_bgc_id:
        raise ValueError("Only match was to self")
    if len(ids) == 0:
        raise ValueError("No matches found")

    assemblies = [search_object.input_assembly] + list(assemblies)[:50]
    z = SerializeToClustermap(
        sg_object=search_object.sg_object,
        sorted_bgcs=assemblies,
        link_df=search_object.link_df,
        group_df=search_object.group_df,
    )
    z.write(outpath_clinker)
    return search_object
