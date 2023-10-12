#!/usr/bin/env python3


import argparse
import logging
import os

from socialgene.utils.run_subprocess import run_subprocess

# mmseqs must be present on PATH


parser = argparse.ArgumentParser(
    description="Run MMseqs2 createdb on an input fasta file"
)

parser.add_argument(
    "input_fasta_path",
    help="Path to fasta file",
)

parser.add_argument(
    "target_database_path",
    help="Path to fasta file",
)


def create_database(fasta_path, database_path):
    """Run MMseqs2 createdb on an input fasta file

    Args:
        fasta_path (str): Path to fasta file
        target_database (str): Path to MMseqs2 database
    """

    command_list = [
        "mmseqs",
        "createdb",
        str(fasta_path),
        str(database_path),
    ]
    mes = run_subprocess(
        command_list=command_list, check=False, shell=False, capture_output=True
    )
    logging.debug(mes)
    if not os.path.exists(database_path):
        logging.warning("No output from MMseqs2 create_database")


def main():
    args = parser.parse_args()

    fasta_path = args.input_fasta_path
    target_database = args.target_database_path

    if not os.path.exists(fasta_path):
        raise FileExistsError(fasta_path)

    create_database(fasta_path, target_database)


if __name__ == "__main__":
    main()
