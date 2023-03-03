# python dependencies
import argparse
import glob
from pathlib import Path
import csv

# external dependencies

# internal dependencies
from socialgene.base.socialgene import SocialGene
from socialgene.utils.logging import log
from socialgene.config import env_vars
from socialgene.utils.writers import write_tsv


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
    "--collect_tables_in_memory",
    metavar="bool",
    help="Should all tables be collected in RAM before writing to disk?",
    required=False,
    default=False,
    action=argparse.BooleanOptionalAction,
)


def fasta_gen(obj):
    for k, v in obj.proteins.items():
        yield f">{k}\n{v.sequence}"


def export_tables(
    sequence_files_glob, outdir, n_fasta_splits, collect_tables_in_memory
):
    outdir = Path(outdir)
    # find files using the input glob
    sequence_files = glob.glob(sequence_files_glob)
    sequence_files = (Path(i).resolve() for i in sequence_files)
    if not sequence_files:
        raise IOError("No matching file(s)")

    if collect_tables_in_memory:
        # Avoid overhead of opening file connections per each file
        # Instead accumulate tables in memory and dump at end
        # still parse into individual SocialGene() object to avoid overhead of
        # appending to large lists
        print("hi1")
        fasta_accum = []
        protein_info_table = []
        assembly_to_locus_table = []
        loci_table = []
        assembly_table = []
        assembly_to_taxid_table = []

        for i in sequence_files:
            socialgene_object = SocialGene()
            socialgene_object.parse(Path(i))
            fasta_accum.append(fasta_gen(socialgene_object))
            protein_info_table.append(socialgene_object.protein_info_table())
            assembly_to_locus_table.append(socialgene_object.assembly_to_locus_table())
            loci_table.append(socialgene_object.loci_table())
            assembly_table.append(socialgene_object.assembly_table())
            assembly_to_taxid_table.append(socialgene_object.assembly_to_taxid_table())
        with open("fasta.faa", "w") as handle:
            for i in fasta_accum:
                handle.writelines("{}\n".format(x) for x in i)
        for i in SocialGene.tsv_tablenames():
            with open(i.removesuffix("_table"), "w") as handle:
                for i in locals()[i]:
                    tsv_writer = csv.writer(
                        handle, delimiter="\t", quotechar='"', quoting=csv.QUOTE_MINIMAL
                    )
                    for ii in i:
                        tsv_writer.writerow(ii)
    else:
        print("hi2")
        # Export tables after parsing each file
        for i in sequence_files:
            print("hi3")
            socialgene_object = SocialGene()
            socialgene_object.parse(Path(i))
            socialgene_object.write_n_fasta(
                outdir=outdir, n_splits=n_fasta_splits, mode="a"
            )
            for i in SocialGene.tsv_tablenames():
                socialgene_object.write_table(
                    outdir=outdir, type=i, filename=i, mode="a"
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
