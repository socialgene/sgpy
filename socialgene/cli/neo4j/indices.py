import argparse
from socialgene.neo4j.schema.graph_schema import GraphSchema
from socialgene.utils.logging import log

from socialgene.dbmodifiers.massage.indices import (
    assembly_uid,
    hmm_uid,
    nucleotide_external_id,
    nucleotide_uid,
    protein_uid,
    taxonomy_uid,
)

parser = argparse.ArgumentParser(
    description="Add indices to a SocialGene Neo4j Database"
)
parser.add_argument(
    "--all",
    help="",
    default=False,
    required=False,
    action="store_true",
)
parser.add_argument(
    "--labels",
    help="Node labels to try and add indices to",
    default=False,
    required=False,
    nargs="+",
)

def main():
    args = parser.parse_args()
    as_dict={i.neo4j_label:i for i in GraphSchema.ALL_NODES}
    if not any([i for i in args.__dict__.values()]):
        parser.print_help()
        print(f"Available labels: {sorted(list(as_dict.keys()))}")
        return
    if args.all:
        for i in GraphSchema.ALL_NODES:
            i().add_constraints_to_neo4j()
    elif args.labels:
        for i in args.labels:
            if i in as_dict:
                as_dict[i]().add_constraints_to_neo4j()
            else:
                log.warning(f"Label {i} not found in GraphSchema.ALL_NODES")
                log.warning(f"Available labels: {as_dict.keys()}")

if __name__ == "__main__":
    main()
