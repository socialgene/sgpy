import argparse

from socialgene.neo4j.schema.socialgene_modules import SocialgeneModules
from socialgene.utils.logging import log

parser = argparse.ArgumentParser(
    description="Export header files for Neo4j admin import"
)

parser.add_argument(
    "--outdir",
    help="Output directory filepath",
    required=True,
)
parser.add_argument(
    "--sg_modules",
    help="Modules for data import; current choices: ['base', 'blastp', 'hmms', 'mmseqs2', 'ncbi_taxonomy', 'paired_omics']",
    required=True,
    nargs="+",
)


def main():
    args = parser.parse_args()
    log.info(args.sg_modules)
    sg_mod = SocialgeneModules()
    sg_mod.add_modules(args.sg_modules)
    sg_mod.write_neo4j_headers(outdir=args.outdir)


if __name__ == "__main__":
    main()
