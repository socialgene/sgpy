import argparse

from socialgene.dbmodifiers.massage.massage import (
    add_antismash_regions_as_edges,
    add_antismash_regions_as_nodes,
    add_protein_descriptions,
    add_taxonomic_name_to_assembly,
    culture_collections_as_nodes_rels,
    fix_mibig_taxonomy,
    set_mibig_bgc,
)
from socialgene.utils.logging import log

NCBI_RANKS = [
    "biotype",
    "clade",
    "class",
    "cohort",
    "family",
    "forma",
    "forma specialis",
    "genotype",
    "genus",
    "infraclass",
    "infraorder",
    "isolate",
    "kingdom",
    "morph",
    "no rank",
    "order",
    "parvorder",
    "pathogroup",
    "phylum",
    "section",
    "series",
    "serogroup",
    "serotype",
    "species",
    "species group",
    "species subgroup",
    "strain",
    "subclass",
    "subcohort",
    "subfamily",
    "subgenus",
    "subkingdom",
    "suborder",
    "subphylum",
    "subsection",
    "subspecies",
    "subtribe",
    "superclass",
    "superfamily",
    "superkingdom",
    "superorder",
    "superphylum",
    "tribe",
    "varietas",
]
parser = argparse.ArgumentParser(description="Parse NcbiAssembliessdsd taxonomy")

parser.add_argument(
    "--antismash_as_edges",
    help="",
    default=False,
    required=False,
    action=argparse.BooleanOptionalAction,
)
parser.add_argument(
    "--antismash_as_nodes",
    help="",
    default=False,
    required=False,
    action=argparse.BooleanOptionalAction,
)
parser.add_argument(
    "--protein_descriptions",
    help="",
    default=False,
    required=False,
    action=argparse.BooleanOptionalAction,
)
parser.add_argument(
    "--fix_mibig_taxonomy",
    help="",
    default=False,
    required=False,
    action=argparse.BooleanOptionalAction,
)
parser.add_argument(
    "--set_mibig_bgc",
    help="",
    default=False,
    required=False,
    action=argparse.BooleanOptionalAction,
)
parser.add_argument(
    "--culture_collections",
    help="",
    default=False,
    required=False,
    action=argparse.BooleanOptionalAction,
)

parser.add_argument(
    "--taxa",
    nargs="+",
    metavar="N",
    type=str,
    help="Adds taxon as a node property of the given rank for each assembly (eg: genus family)",
    required=False,
    default=None,
    choices=NCBI_RANKS,
)


def main():
    args = parser.parse_args()
    # defaults are all negative, print help if nothing is provided
    if not any([i for i in args.__dict__.values()]):
        parser.print_help()
        return
    if args.antismash_as_edges:
        add_antismash_regions_as_edges()
    if args.antismash_as_nodes:
        add_antismash_regions_as_nodes()
    if args.protein_descriptions:
        add_protein_descriptions()
    if args.fix_mibig_taxonomy:
        fix_mibig_taxonomy()
    if args.set_mibig_bgc:
        set_mibig_bgc()
    if args.culture_collections:
        culture_collections_as_nodes_rels()
    if args.taxa:
        for i in args.taxa:
            if i in NCBI_RANKS:
                add_taxonomic_name_to_assembly(i)
            else:
                log.warning(f"{i} not in the list of known ranks")


if __name__ == "__main__":
    main()
