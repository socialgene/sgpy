from typing import List

import argparse
import glob
from pathlib import Path

from socialgene.base.socialgene import SocialGene
from socialgene.config import env_vars
from socialgene.utils.logging import log

parser = argparse.ArgumentParser(
    description="Export Socialagene object into TSV files for Nextflow/Neo4j"
)
parser.add_argument(
    "--sequence_files_glob",
    metavar="filepath",
    help='A "glob" style input path. See https://docs.python.org/3/library/glob.html for info about globs',
    required=False,
)
parser.add_argument(
    "-f",
    "--sequence_file",
    help="Filepath, may be repeated to add multiple paths",
    nargs="+",
    action="append",
    required=False,
)
parser.add_argument(
    "--outdir",
    metavar="filepath",
    help="Output directory filepath",
    required=True,
)
parser.add_argument(
    "--n_fasta_splits",
    metavar="filepath",
    help="Number of fasta files to create",
    required=True,
    default=1,
)
parser.add_argument(
    "--collect_tables_in_memory",
    metavar="bool",
    help="Should all tables be collected in RAM before writing to disk?",
    required=False,
    default=False,
    action=argparse.BooleanOptionalAction,
)


def _writer(socialgene_object, outdir, n_fasta_splits):
    socialgene_object.write_n_fasta(outdir=outdir, n_splits=n_fasta_splits, mode="a")
    for i in SocialGene._genomic_info_export_tablenames:
        socialgene_object.write_table(
            outdir=outdir, tablename=i, filename=i.removesuffix("_table"), mode="a"
        )


def export_tables(
    outdir: str,
    n_fasta_splits: int,
    collect_tables_in_memory: bool,
    file_list: List = None,
    sequence_files_glob: str = "",
):
    outdir = Path(outdir)
    sequence_files = list()
    # find files using the input glob
    if sequence_files_glob:
        sequence_files.extend(glob.glob(sequence_files_glob))
    # and/or use input file paths
    if file_list:
        sequence_files.extend(file_list)
    sequence_files = (Path(i).resolve() for i in sequence_files)
    if not sequence_files:
        raise IOError("No matching file(s)")
    # Import all files into the same SG object
    socialgene_object = SocialGene()
    for i in sequence_files:
        socialgene_object.parse(Path(i))
        if not collect_tables_in_memory:
            _writer(
                socialgene_object=socialgene_object,
                outdir=outdir,
                n_fasta_splits=n_fasta_splits,
            )
            socialgene_object = SocialGene()
    if collect_tables_in_memory:
        _writer(
            socialgene_object=socialgene_object,
            outdir=outdir,
            n_fasta_splits=n_fasta_splits,
        )


def main():
    args = parser.parse_args()
    log.info(f"Socialgene variables: \n{env_vars}")
    export_tables(
        sequence_files_glob=args.sequence_files_glob,
        outdir=args.outdir,
        n_fasta_splits=int(args.n_fasta_splits),
        collect_tables_in_memory=args.collect_tables_in_memory,
    )


if __name__ == "__main__":
    main()
