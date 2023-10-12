import os
from pathlib import Path

from socialgene.cli.nextflow.export_protein_loci_assembly_tables import export_tables
from socialgene.config import env_vars

FIXTURE_DIR = os.path.dirname(os.path.realpath(__file__))
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)

TEST_GENOMES_DIR = os.path.join(FIXTURE_DIR, "python", "data", "test_genomes")


def doit(EXPECTED_RESULTS_DIR):
    if not os.path.exists(EXPECTED_RESULTS_DIR):
        os.makedirs(EXPECTED_RESULTS_DIR)
    genome_filepaths = [
        Path(TEST_GENOMES_DIR, i)
        for i in [
            "BGC0000493.gbk",
            "BGC0000492.gbk",
            "lagriamide_mibig_bgc0001946.gbk",
            "GCA_028106715.1_ASM2810671v1_genomic.gbff.gz",
        ]
    ]
    export_tables(
        file_list=genome_filepaths,
        outdir=Path(EXPECTED_RESULTS_DIR),
        n_fasta_splits=1,
        collect_tables_in_memory=True,
    )


def main():
    env_vars["HASHING_ALGORITHM"] = "crc64"
    EXPECTED_RESULTS_DIR = os.path.join(
        FIXTURE_DIR, "python", "data", "truths", "parse_using_crc64"
    )
    doit(EXPECTED_RESULTS_DIR)

    env_vars["HASHING_ALGORITHM"] = "crc32"
    EXPECTED_RESULTS_DIR = os.path.join(
        FIXTURE_DIR, "python", "data", "truths", "crc32"
    )
    doit(EXPECTED_RESULTS_DIR)


if __name__ == "__main__":
    main()
