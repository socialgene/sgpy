import argparse
from pathlib import Path
from typing import List

from socialgene.clustermap.serialize import SerializeToClustermap
from socialgene.config import env_vars
from socialgene.search.hmmer import SearchDomains
from socialgene.utils.logging import log

env_vars["NEO4J_URI"] = "bolt://localhost:7687"

parser = argparse.ArgumentParser(
    description="Search a SocialGene database for input gene clusters similar to an input gene cluster"
)


parser.add_argument(
    "--gbk_path",
    metavar="input_filepath",
    type=Path,
    help="Path to the query gene cluster GenBank-format file",
    required=True,
)
parser.add_argument(
    "--hmm_dir",
    metavar="hmm_directory",
    type=str,
    help="Path to the directory containing HMM files created by SocialGene's Nextflow workflow",
    required=True,
)
parser.add_argument(
    "--assemblies_must_have_x_matches",
    metavar="fraction",
    type=float,
    default=0.6,
    help="Minimum query proteins an assembly must have (<1 == fraction of query proteins)",
)
parser.add_argument(
    "--nucleotide_sequences_must_have_x_matches",
    metavar="fraction",
    type=float,
    default=0.6,
    help="Minimum query proteins a nucleotide sequence must have (<1 == fraction of query proteins)",
)
parser.add_argument(
    "--gene_clusters_must_have_x_matches",
    metavar="fraction",
    type=float,
    default=0.6,
    help="Minimum query proteins a target gene cluster sequence must have (<1 == fraction of query proteins)",
)
parser.add_argument(
    "--break_bgc_on_gap_of",
    metavar="value",
    type=int,
    default=10000,
    help="Split gene clusters if gap between matched proteins is greater than this value",
)
parser.add_argument(
    "--target_bgc_padding",
    metavar="value",
    type=int,
    default=2000,
    help="Pull proteins x-nucleotides on either side of target gene cluster",
)
parser.add_argument(
    "--max_domains_per_protein",
    metavar="value",
    type=int,
    default=3,
    help="Maximum number of domains per protein used in search; 0 for all",
)
parser.add_argument(
    "--max_outdegree",
    metavar="value",
    type=int,
    help="Maximum outdegree a query domain can have; 0 for all",
)
parser.add_argument(
    "--max_query_proteins",
    metavar="value",
    type=int,
    default=10,
    help="Maximum number of query proteins to use in search; 0 for all",
)
parser.add_argument(
    "--scatter",
    metavar="value",
    help="Pull proteins from across the width of the input gene cluster",
    action=argparse.BooleanOptionalAction,
    default=False,
)
parser.add_argument(
    "--locus_tag_bypass_list",
    metavar="value",
    nargs="+",
    help="List of locus tags to bypass search filter",
)
parser.add_argument(
    "--protein_id_bypass_list",
    metavar="value",
    nargs="+",
    help="List of protein IDs to bypass search filter",
)
parser.add_argument(
    "--only_culture_collection",
    metavar="value",
    help="Only search genomes from culture collections",
    action=argparse.BooleanOptionalAction,
    default=False,
)
parser.add_argument(
    "--frac",
    metavar="fraction",
    type=float,
    default=0.75,
    help="Fraction of domains that equals a match in the initial database scan",
)
parser.add_argument(
    "--run_async",
    metavar="value",
    help="Run the initial database scan using asynchronously",
    action=argparse.BooleanOptionalAction,
    default=True,
)
parser.add_argument(
    "--analyze_with",
    metavar="value",
    type=str,
    default="hmmer",
    help="Tool to use for reciprocal best hit comparison; one of 'hmmer', 'blastp', 'mmseqs2'",
)
args = parser.parse_args()


def search_bgc(
    gbk_path: str,
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
):
    log.info(f"Running search with args: {locals()}")
    search_object = SearchDomains(
        gbk_path=gbk_path,
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
        tool=analyze_with, argstring="--fast --max-hsps 1", cpus=1
    )
    search_object._choose_group()
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
    z.write("/Users/chase/Documents/test/cmap/clinker/clinker/plot/data.json")


def main():
    args = parser.parse_args()

    if args.max_domains_per_protein == 0:
        args.max_domains_per_protein = None
    if args.max_outdegree == 0:
        args.max_outdegree = None
    if args.max_query_proteins == 0:
        args.max_query_proteins = None
    _ = search_bgc(
        gbk_path=args.gbk_path,
        hmm_dir=args.hmm_dir,
        use_neo4j_precalc=True,
        assemblies_must_have_x_matches=args.assemblies_must_have_x_matches,
        nucleotide_sequences_must_have_x_matches=args.nucleotide_sequences_must_have_x_matches,
        gene_clusters_must_have_x_matches=args.gene_clusters_must_have_x_matches,
        break_bgc_on_gap_of=args.break_bgc_on_gap_of,
        target_bgc_padding=args.target_bgc_padding,
        max_domains_per_protein=args.max_domains_per_protein,
        max_outdegree=args.max_outdegree,
        max_query_proteins=args.max_query_proteins,
        scatter=args.scatter,
        locus_tag_bypass_list=args.locus_tag_bypass_list,
        protein_id_bypass_list=args.protein_id_bypass_list,
        only_culture_collection=args.only_culture_collection,
        frac=args.frac,
        run_async=args.run_async,
    )


if __name__ == "__main__":
    main()
