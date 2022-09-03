#!/usr/bin/env python

# python dependencies
import argparse
import glob
from pathlib import Path
from multiprocessing import Pool, cpu_count

# external dependencies
from rich.progress import Progress

# internal dependencies
from socialgene.base.socialgene import SocialGene
from socialgene.utils.logging import log
from socialgene.config import env_vars
from socialgene.utils.chunker import chunk_a_list_with_numpy


parser = argparse.ArgumentParser(
    description="Export Socialagene object into TSV files for Nextflow/Neo4j"
)
parser.add_argument(
    "--sequence_files_glob",
    metavar="filepath",
    help='A "glob" style input path. See https://docs.python.org/3/library/glob.html for info about globs',
    required=True,
    # nargs="+", not used since moved to glob
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
    "--parallel",
    metavar="bool",
    help="",
    required=False,
    default=False,
    action=argparse.BooleanOptionalAction,
)


def simple_parse(input_pathlist):
    sg_object = SocialGene()
    for input_path in input_pathlist:
        try:
            sg_object.parse(input_path, keep_sequence=False)
        except Exception as e:
            log.error(e)
    return sg_object


def export_tables(sequence_files_glob, outdir, n_fasta_splits, parallel=False, cpus=0):
    outdir = Path(outdir)
    socialgene_object = SocialGene()
    # find files using the input glob
    sequence_files = glob.glob(sequence_files_glob)
    sequence_files = [Path(i).resolve() for i in sequence_files]
    if sequence_files == []:
        raise IOError("No matching file(s)")
    if isinstance(sequence_files, list):
        if parallel:
            if not cpus == 0:
                cpus = cpus
            elif cpu_count() < 2:
                cpus = 1
            else:
                cpus = cpu_count() - 1
            chunked_list = chunk_a_list_with_numpy(
                input_list=sequence_files, n_chunks=cpus
            )
            chunked_list = [i for i in chunked_list]
            with Progress(transient=True) as progress:
                task = progress.add_task("Reading...", total=len(chunked_list))
                with Pool(processes=cpus) as p:
                    imap_unordered_it = p.imap_unordered(simple_parse, chunked_list)
                    for i in imap_unordered_it:
                        i.export_locus_to_protein(outdir=outdir)
                        i.export_protein_info(outdir=outdir)
                        i.export_loci(outdir=outdir)
                        i.export_assembly(outdir=outdir)
                        i.export_assembly_to_locus(outdir=outdir)
                        i.export_assemby_to_taxid(outdir=outdir)
                        progress.update(task, advance=1)

        else:
            for i in sequence_files:
                log.info(i)
                socialgene_object.parse(Path(i))
    else:
        raise TypeError("Unexpected sequence_files type")
    socialgene_object.write_n_fasta(outdir=outdir, n_splits=n_fasta_splits)
    socialgene_object.export_locus_to_protein(outdir=outdir)
    socialgene_object.export_protein_info(outdir=outdir)
    socialgene_object.export_loci(outdir=outdir)
    socialgene_object.export_assembly(outdir=outdir)
    socialgene_object.export_assembly_to_locus(outdir=outdir)
    socialgene_object.export_assemby_to_taxid(outdir=outdir)
    # socialgene_object.write_n_fasta(outdir=outdir, n_splits=n_fasta_splits)


def main():
    args = parser.parse_args()
    log.info(f"Socialgene variables: \n{env_vars}")
    export_tables(
        sequence_files_glob=args.sequence_files_glob,
        outdir=args.outdir,
        n_fasta_splits=int(args.n_fasta_splits),
        parallel=args.parallel,
    )


if __name__ == "__main__":
    main()
