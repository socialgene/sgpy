#!/usr/bin/env python3


import argparse
import logging
import os

from socialgene.utils.run_subprocess import run_subprocess

# mmseqs must be present on PATH


parser = argparse.ArgumentParser(
    description="Run MMseqs2 subsetdb on an input fasta file"
)

parser.add_argument(
    "idfile",
    help="Path to idfile file",
)

parser.add_argument(
    "olddb",
    help="Path to olddb",
)

parser.add_argument(
    "newdb",
    help="Path to newdb",
)


def createsubdb(olddb, newdb, idfile):
    command_list = [
        "mmseqs",
        "createsubdb",
        str(idfile),
        str(olddb),
        str(newdb),
    ]
    mes = run_subprocess(
        command_list=command_list, check=False, shell=False, capture_output=True
    )
    logging.debug(mes)
    if not os.path.exists(newdb):
        logging.warning("No output from MMseqs2 createsubdb")


def main():
    args = parser.parse_args()
    createsubdb(args.olddb, args.newdb, args.idfile)


if __name__ == "__main__":
    main()
