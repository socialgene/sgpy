import os
import tempfile
from pathlib import Path

import pandas as pd

from socialgene.base.molbio import Protein
from socialgene.compare_proteins.base import BlastTab, BlastTab_COLUMNS
from socialgene.utils.logging import log
from socialgene.utils.run_subprocess import run_subprocess


class DiamondBlastp(BlastTab):
    def __init__(self):
        super().__init__()
        self.temp_db_path = None
        self.df = None
        self.name = "Diamond BLASTp"

    def make_db(
        self,
        fasta_string: str,
        db_path=None,
        threads=1,
        **kwargs,
    ):
        db_path = Path(db_path).with_suffix(".dmnd")
        # run diamond makedb on fasta_file allow passing kwargs to diamond
        command_list = [
            "diamond",
            "makedb",
            "--db",
            str(db_path),
            "--threads",
            str(threads),
            "--verbose",
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
                f"No database from DiamondBlastp make_db found at {db_path}"
            )

    @staticmethod
    def run(fasta_path, db_path=None, threads=1, argstring="", **kwargs):
        # run diamond blastp on query_fasta_string against db_path
        if not Path(db_path).exists():
            raise FileNotFoundError(
                f"No database from Diamond makedb found at {db_path}"
            )
        with tempfile.TemporaryDirectory() as tmpdirname:
            outpath = os.path.join(tmpdirname, "result.m8")
            command_list = [
                "diamond",
                "blastp",
                "--query",
                str(fasta_path),
                "--db",
                str(db_path),
                "--out",
                str(outpath),
                "--threads",
                str(threads),
                "--verbose",
                "--outfmt",
                "6",
                "qseqid",
                "sseqid",
                "pident",
                "length",
                "mismatch",
                "gapopen",
                "qstart",
                "qend",
                "sstart",
                "send",
                "evalue",
                "bitscore",
                "qlen",
                "slen",
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
            log.debug(mes)
            if not os.path.exists(outpath):
                log.warning("No output from Diamond search")
            else:
                # sorted for reproducibility
                return (
                    pd.read_csv(
                        outpath,
                        sep="\t",
                        header=None,
                        names=BlastTab_COLUMNS,
                        dtype=BlastTab_COLUMNS,
                    )
                    .sort_values(["query", "score"], ascending=False)
                    .reset_index(drop=True)
                )

    def compare_proteins(
        self, queries, targets, cpus=1, argstring="--fast --no-self-hits --max-hsps 1"
    ):
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
            target_dmnd_path = Path(tmpdirname, "targets.dmnd")
            with open(query_fasta_path, "w") as h:
                for protein in queries:
                    h.write(protein.fasta_string_defline_uid)
            self.make_db(
                fasta_string="".join([i.fasta_string_defline_uid for i in targets]),
                db_path=target_dmnd_path,
            )
            return self.run(
                fasta_path=query_fasta_path,
                db_path=target_dmnd_path,
                threads=cpus,
                argstring=argstring,
            )
