import argparse
import glob
from pathlib import Path

import socialgene.utils.protein_sqlite as ps
from socialgene.utils.logging import log

parser = argparse.ArgumentParser(
    description="Create a sqlite database from multiple tsv files (hashid\taccession)"
)
parser.add_argument(
    "--input",
    metavar="filepath",
    help='A "glob" style input path. Make sure to use single quotes (\'), not double("). See https://docs.python.org/3/library/glob.html for info about globs.',
    required=True,
)
parser.add_argument(
    "--outpath",
    metavar="filepath",
    help="Path sqlite file will be written to",
    required=True,
)
parser.add_argument(
    "--build_in_memory",
    metavar="bool",
    help="Path sqlite file will be written to",
    required=False,
    default=True,
    action=argparse.BooleanOptionalAction,
)


def main():
    args = parser.parse_args()

    paths = glob.glob(args.input)
    paths = [Path(i) for i in paths]
    log.info(f"Found {len(paths)} files to insert into {args.outpath}")
    log.info("Note: database is built")
    log.info(args.outpath)
    ps.create_db(tsv_list=paths, outpath=args.outpath, inmem=args.build_in_memory)


if __name__ == "__main__":
    main()
