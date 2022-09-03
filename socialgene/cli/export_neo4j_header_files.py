#!/usr/bin/env python

# python dependencies
import argparse

# external dependencies

# internal dependencies
from socialgene.neo4j.sg_modules import write_neo4j_headers, parse_hmmlist_input
from socialgene.utils.logging import log

parser = argparse.ArgumentParser(
    description="Export header files for Neo4j admin import"
)

parser.add_argument(
    "--outdir",
    metavar="filepath",
    help="Output directory filepath",
    required=True,
)
parser.add_argument(
    "--sg_modules",
    metavar="list",
    help="Modules for data import; current choices: ['base', 'blastp', 'hmms', 'mmseqs2', 'ncbi_taxonomy', 'paired_omics']",
    required=True,
    nargs="+",
)

parser.add_argument(
    "--hmmlist",
    metavar="comma sep list",
    help="List of HMM models used",
    required=True,
    default="base",
    nargs="+",
)


def main():
    args = parser.parse_args()
    log.info(args.sg_modules)
    log.info(args.hmmlist)
    parsed_hmmlist = parse_hmmlist_input(args.hmmlist)
    log.info(parsed_hmmlist)
    write_neo4j_headers(
        sg_modules=args.sg_modules, hmmlist=parsed_hmmlist, outdir=args.outdir
    )


if __name__ == "__main__":
    main()
