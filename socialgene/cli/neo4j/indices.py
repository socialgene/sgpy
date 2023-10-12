import argparse

from socialgene.dbmodifiers.massage.indices import (
    assembly_uid,
    hmm_uid,
    nucleotide_external_id,
    nucleotide_uid,
    protein_uid,
)

parser = argparse.ArgumentParser(
    description="Add indices to a SocialGene Neo4j Database"
)
parser.add_argument(
    "--all",
    help="",
    default=False,
    required=False,
    action=argparse.BooleanOptionalAction,
)
parser.add_argument(
    "--assembly_uid",
    help="",
    default=False,
    required=False,
    action=argparse.BooleanOptionalAction,
)
parser.add_argument(
    "--protein_uid",
    help="",
    default=False,
    required=False,
    action=argparse.BooleanOptionalAction,
)

parser.add_argument(
    "--hmm_uid",
    help="",
    default=False,
    required=False,
    action=argparse.BooleanOptionalAction,
)
parser.add_argument(
    "--nucleotide_uid",
    help="",
    default=False,
    required=False,
    action=argparse.BooleanOptionalAction,
)
parser.add_argument(
    "--nucleotide_external_id",
    help="",
    default=False,
    required=False,
    action=argparse.BooleanOptionalAction,
)


def main():
    args = parser.parse_args()
    # defaults are all negative, print help if nothing is provided
    if not any([i for i in args.__dict__.values()]):
        parser.print_help()
        return
    if args.assembly_uid or args.all:
        assembly_uid()
    if args.protein_uid or args.all:
        protein_uid()
    if args.hmm_uid or args.all:
        hmm_uid()
    if args.nucleotide_uid or args.all:
        nucleotide_uid()
    if args.nucleotide_external_id or args.all:
        nucleotide_external_id()


if __name__ == "__main__":
    main()
