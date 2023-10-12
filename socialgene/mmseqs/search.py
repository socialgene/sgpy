#!/usr/bin/env python3


import argparse
import logging
import os
import tempfile

import pandas as pd

from socialgene.utils.run_subprocess import run_subprocess

# mmseqs must be present on PATH


parser = argparse.ArgumentParser(
    description="Run MMseqs2 search on an input amino acid sequence"
)

parser.add_argument(
    "input_fasta_path",
    help="Path to fasta file",
)

parser.add_argument(
    "target_database_path",
    help="Path to fasta file",
)

MMSEQS_OUT_COLUMNS = {
    "query": str,
    "target": str,
    "pident": "Float32",
    "alnlen": "Int32",
    "mismatch": "Int32",
    "gapopen": "Int32",
    "qstart": "Int32",
    "qend": "Int32",
    "tstart": "Int32",
    "tend": "Int32",
    "evalue": "Float32",
    "bits": "Int32",
}


def search(fasta_path, target_database, argstring=""):
    """Search an input fasta file against a target database with external MMseqs2 program

        Args:
            fasta_path (str): Path to fasta file
            target_database (str): Path to MMseqs2 database
            argstring (str): additonal arguments to pass to mmseqs search, as a string
    =
        Returns:
            pandas_df: dataframe of the mmseqs search result
    """
    with tempfile.TemporaryDirectory() as tmpdirname:
        outpath = os.path.join(tmpdirname, "result.m8")
        command_list = [
            "mmseqs",
            "easy-search",
            str(fasta_path),
            str(target_database),
            str(outpath),
            str(tmpdirname),
            "--format-mode",
            "0",
        ]
        if argstring:
            argstring = [str(i) for i in argstring.split()]
            command_list.extend(argstring)
        mes = run_subprocess(
            command_list=command_list, check=False, shell=False, capture_output=True
        )
        logging.info(mes)
        if not os.path.exists(outpath):
            logging.warning("No output from MMseqs2 search")
        else:
            return (
                pd.read_csv(
                    outpath,
                    sep="\t",
                    names=MMSEQS_OUT_COLUMNS,
                    dtype=MMSEQS_OUT_COLUMNS,
                )
                .sort_values(["pident", "bits"], ascending=False)
                .reset_index(drop=True)
            )


def main():
    args = parser.parse_args()

    fasta_path = args.input_fasta_path
    target_database = args.target_database_path

    if not os.path.exists(fasta_path):
        raise FileExistsError(fasta_path)
    if not os.path.exists(target_database):
        raise FileExistsError(target_database)

    print(search(fasta_path, target_database))


if __name__ == "__main__":
    main()
