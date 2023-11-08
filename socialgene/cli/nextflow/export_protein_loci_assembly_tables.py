import argparse
import glob
from pathlib import Path
from typing import List

from socialgene.base.socialgene import SocialGene
from socialgene.config import env_vars
from socialgene.utils.logging import log

parser = argparse.ArgumentParser(
    description="Export Socialagene object into TSV files for Nextflow/Neo4j"
)
parser.add_argument(
    "--sequence_files_glob",
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
    help="Output directory filepath",
    required=True,
)
parser.add_argument(
    "--n_fasta_splits",
    help="Number of fasta files to create",
    required=True,
    default=1,
)
parser.add_argument(
    "--include_sequences",
    help="Should sequences be included in the database?",
    required=False,
    default=False,
    action=argparse.BooleanOptionalAction,
)
parser.add_argument(
    "--collect_tables_in_memory",
    help="Should all tables be collected in RAM before writing to disk?",
    required=False,
    default=False,
    action=argparse.BooleanOptionalAction,
)
parser.add_argument(
    "--compression",
    help="Gzip compress output?",
    required=False,
    default=False,
    action=argparse.BooleanOptionalAction,
)
parser.add_argument(
    "--defline_magic",
    help="Parse out fasta deflines with defined structures, currently only works for uniprot (deflines that begin with: sp| or tr|)",
    type=bool,
    required=False,
    default=False,
)


def _writer(socialgene_object, outdir, n_fasta_splits, compression, include_sequences):
    socialgene_object.write_n_fasta(
        outdir=outdir, n_splits=n_fasta_splits, mode="a", compression=compression
    )
    for i in SocialGene._export_table_names:
        socialgene_object.write_table(
            outdir=outdir,
            tablename=i,
            filename=i.removeprefix("table_"),
            mode="a",
            compression=compression,
            include_sequences=include_sequences,
        )


def export_tables(
    outdir: str,
    n_fasta_splits: int,
    collect_tables_in_memory: bool,
    compression: str,
    file_list: List = None,
    sequence_files_glob: str = "",
    include_sequences: bool = False,
    defline_magic: bool = False,
):
    """
    The function `export_tables` takes in various parameters, including an output directory, the number
    of fasta splits, a flag to determine if tables should be collected in memory, a compression method,
    and a list of files or a glob pattern to find files. It then imports the files into a SocialGene
    object and writes the tables to the specified output directory.

    Args:
      outdir (str): The `outdir` parameter is a string that specifies the directory where the exported
    tables will be saved.
      n_fasta_splits (int): The `n_fasta_splits` parameter specifies the number of splits to be made
    when exporting the tables. It is used by the `_writer` function to determine how many splits to
    create when writing the tables to files.
      collect_tables_in_memory (bool): A boolean flag indicating whether the tables should be collected
    in memory before writing them to files. If set to True, the tables will be stored in memory and then
    written to files. If set to False, the tables will be written to files immediately after parsing
    each input file.
      compression (str): The "compression" parameter is a string that specifies the compression format
    to use when exporting the tables. It can be one of the following options: "gzip", "bz2", "lzma", or
    "none".
      file_list (List): A list of file paths to input sequence files.
      sequence_files_glob (str): A string representing a glob pattern used to find sequence files. This
    pattern can include wildcards (*) to match multiple files.
      include_sequences (bool): A boolean flag indicating whether the sequences should be included in the database
      defleine_magic (bool): Parse out fasta deflines with defined structures, currently only works for uniprot (deflines that begin with: sp| or tr|)

    """
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
    if collect_tables_in_memory:
        # Import all files into the same SG object
        socialgene_object = SocialGene()
        for i in sequence_files:
            socialgene_object.parse(Path(i), defline_magic=defline_magic)
        _writer(
            socialgene_object=socialgene_object,
            outdir=outdir,
            n_fasta_splits=n_fasta_splits,
            compression=compression,
            include_sequences=include_sequences,
        )
    else:
        for i in sequence_files:
            socialgene_object = SocialGene()
            socialgene_object.parse(Path(i), defline_magic=defline_magic)
            _writer(
                socialgene_object=socialgene_object,
                outdir=outdir,
                n_fasta_splits=n_fasta_splits,
                compression=compression,
                include_sequences=include_sequences,
            )


def main():
    args = parser.parse_args()
    log.info(f"Socialgene variables: \n{env_vars}")
    if args.compression:
        compression = "gzip"
    else:
        compression = None
    export_tables(
        sequence_files_glob=args.sequence_files_glob,
        outdir=args.outdir,
        n_fasta_splits=int(args.n_fasta_splits),
        collect_tables_in_memory=args.collect_tables_in_memory,
        compression=compression,
        include_sequences=args.include_sequences,
        defline_magic=args.defline_magic,
    )


if __name__ == "__main__":
    main()
