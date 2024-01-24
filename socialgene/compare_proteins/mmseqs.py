#!/usr/bin/env python3

import logging
import os
import tempfile
from pathlib import Path

import pandas as pd

from socialgene.base.molbio import Protein
from socialgene.compare_proteins.base import BlastTab, BlastTab_COLUMNS
from socialgene.utils.logging import log
from socialgene.utils.run_subprocess import run_subprocess

# mmseqs must be present on PATH


class MMseqsEasySearch(BlastTab):
    def __init__(self):
        super().__init__()
        self.temp_db_path = None
        self.df = None
        self.name = "MMseqs2 EasySearch"

    @staticmethod
    def make_db(
        fasta_string: str,
        db_path=None,
        **kwargs,
    ):
        """Run MMseqs2 createdb on an input fasta file

        Args:
            fasta_path (str): Path to fasta file
            target_database (str): Path to MMseqs2 database
        """

        command_list = [
            "mmseqs",
            "createdb",
            "stdin",
            str(db_path),
        ]
        mes = run_subprocess(
            command_list=command_list,
            input=fasta_string,
            text=True,
            capture_output=True,
            status=False,
            **kwargs,
        )
        log.debug(mes)
        if not Path(db_path).exists():
            raise FileNotFoundError(
                f"No database from mmseqs make_db found at {db_path}"
            )

    @staticmethod
    def index_database(database_path):
        """Run MMseqs2 createindex on an MMseqs2 database

        Args:
            target_database (str): Path to MMseqs2 database
        """
        with tempfile.TemporaryDirectory() as tmpdirname:
            command_list = [
                "mmseqs",
                "createindex",
                str(database_path),
                str(tmpdirname),
            ]
            mes = run_subprocess(
                command_list=command_list,
                check=False,
                shell=False,
                capture_output=True,
                status=False,
            )
        logging.debug(mes)
        if not os.path.exists(f"{database_path}.idx"):
            logging.warning("No output from MMseqs2 createindex")

    def run(
        self,
        fasta_path,
        target_database,
        cpus=1,
        argstring="",
    ):
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
                "--format-output",
                "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen",
                "--threads",
                str(cpus),
            ]
            if argstring:
                argstring = [str(i) for i in argstring.split()]
                command_list.extend(argstring)
            mes = run_subprocess(
                command_list=command_list,
                check=False,
                shell=False,
                capture_output=True,
                status=False,
            )
            logging.debug(mes)
            if not os.path.exists(outpath):
                logging.warning("No output from MMseqs2 search")
            else:
                # sorted for reproducibility
                return (
                    pd.read_csv(
                        outpath,
                        sep="\t",
                        names=BlastTab_COLUMNS,
                        dtype=BlastTab_COLUMNS,
                    )
                    .sort_values(
                        [
                            "query",
                            "score",
                        ],
                        ascending=False,
                    )
                    .reset_index(drop=True)
                )

    def compare_proteins(self, queries, targets, cpus=1, argstring="", index=False):
        # loop through protein list and write to temporary fasta file
        if isinstance(queries, Protein):
            queries = [queries]
        if isinstance(targets, Protein):
            targets = [targets]
        if not all([isinstance(i, Protein) for i in queries]):
            raise TypeError("queries must be a list of Protein objects")
        if not all([isinstance(i, Protein) for i in targets]):
            raise TypeError("targets must be a list of Protein objects")

        with tempfile.TemporaryDirectory() as tmpdirname:
            query_fasta_path = Path(tmpdirname, "queries.faa")
            target_dmnd_path = Path(tmpdirname, "targets")
            with open(query_fasta_path, "w") as h:
                for protein in queries:
                    h.write(protein.fasta_string_defline_uid)
            self.make_db(
                fasta_string="".join([i.fasta_string_defline_uid for i in targets]),
                db_path=target_dmnd_path,
            )
            if index:
                self.index_database(target_dmnd_path)
            return self.run(
                fasta_path=query_fasta_path,
                target_database=target_dmnd_path,
                cpus=cpus,
                argstring=argstring,
            )
