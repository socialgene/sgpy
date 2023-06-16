#!/usr/bin/env python3


import argparse
import logging
import os
import tempfile

from socialgene.utils.run_subprocess import run_subprocess

# mmseqs must be present on PATH


parser = argparse.ArgumentParser(
    description="Run MMseqs2 createindex on an MMseqs2 database"
)

parser.add_argument(
    "target_database_path",
    help="Path to fasta file",
)


def index_database(database_path):
    """Run MMseqs2 createindex on an MMseqs2 database

    Args:
        target_database (str): Path to MMseqs2 database
    """
    with tempfile.TemporaryDirectory() as tmpdirname:
        command_list = ["mmseqs", "createindex", database_path, tmpdirname]
        mes = run_subprocess(
            command_list=command_list, check=False, shell=False, capture_output=True
        )
    logging.debug(mes)
    if not os.path.exists(f"{database_path}.idx"):
        logging.warning("No output from MMseqs2 createindex")


def main():
    args = parser.parse_args()
    target_database = args.target_database_path

    if not os.path.exists(target_database):
        raise FileExistsError(target_database)

    index_database(target_database)


if __name__ == "__main__":
    main()
