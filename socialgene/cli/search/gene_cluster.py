import argparse
from pathlib import Path

from socialgene.cli.search.sea import search_bgc
from socialgene.config import env_vars

env_vars["NEO4J_URI"] = "bolt://localhost:7687"

parser = argparse.ArgumentParser(
    description="Search a SocialGene database for gene clusters similar to an input gene cluster"
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
    "--clinker_outpath",
    metavar="output_filepath",
    type=Path,
    help="Path to the output Clinker JSON file",
)
parser.add_argument(
    "--assemblies_must_have_x_matches",
    metavar="fraction",
    type=float,
    default=0.8,
    help="Minimum query proteins an assembly must have (<1 == fraction of query proteins)",
)
parser.add_argument(
    "--nucleotide_sequences_must_have_x_matches",
    metavar="fraction",
    type=float,
    default=0.8,
    help="Minimum query proteins a nucleotide sequence must have (<1 == fraction of query proteins)",
)
parser.add_argument(
    "--gene_clusters_must_have_x_matches",
    metavar="fraction",
    type=float,
    default=0.8,
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
    type=float,
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
    type=float,
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


def main():
    args = parser.parse_args()

    if args.max_domains_per_protein == 0:
        args.max_domains_per_protein = None
    if args.max_outdegree == 0:
        args.max_outdegree = None
    if args.max_query_proteins == 0:
        args.max_query_proteins = None
    _ = search_bgc(
        input=args.gbk_path,
        hmm_dir=args.hmm_dir,
        outpath_clinker=args.clinker_outpath,
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
        analyze_with=args.analyze_with,
    )


if __name__ == "__main__":
    main()
