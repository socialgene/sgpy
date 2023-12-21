import os
import tempfile
from pathlib import Path

import pytest

from socialgene.cli.nextflow.export_protein_loci_assembly_tables import export_tables
from socialgene.config import env_vars

FIXTURE_DIR = os.path.dirname(os.path.realpath(__file__))
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)
TEST_GENOMES_DIR = os.path.join(FIXTURE_DIR, "data", "test_genomes")

FIXTURE_DIR = os.path.dirname(os.path.realpath(__file__))
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)
TRUTHS_DIR = os.path.join(FIXTURE_DIR, "data", "truths", "neo4j_tables")


genome_filepaths = [
    Path(TEST_GENOMES_DIR, i)
    for i in [
        "BGC0000493.gbk",
        "BGC0000492.gbk",
        "lagriamide_mibig_bgc0001946.gbk",
        "GCA_028106715.1_ASM2810671v1_genomic.gbff.gz",
    ]
]


files = [
    "assembly",
    "assembly_to_locus",
    "assembly_to_taxid",
    "fasta_split_1.faa",
    "loci",
    "locus_to_protein",
    "protein_ids",
]
hash_algo = ["crc64", "sha512t24u"]
include_sequences = ["True", "False"]


@pytest.mark.parametrize("item", files)
@pytest.mark.parametrize("hash_algo", hash_algo)
@pytest.mark.parametrize("include_sequence", include_sequences)
def test_gbk_parsing_and_flatfile_for_neo4j_creation_collect_tables_in_memory_false(
    item, hash_algo, include_sequence
):
    env_vars["HASHING_ALGORITHM"] = hash_algo
    with tempfile.TemporaryDirectory() as tmpdirname:
        export_tables(
            file_list=genome_filepaths,
            outdir=tmpdirname,
            n_fasta_splits=1,
            collect_tables_in_memory=True,
            compression=None,
            include_sequences=include_sequence,
        )

        assert Path(tmpdirname, item).exists()
        with open(Path(tmpdirname, item), "r") as h:
            z1 = h.readlines()
            with open(
                Path(TRUTHS_DIR, hash_algo, f"seq_included_{include_sequence}", item),
                "r",
            ) as ht:
                z2 = ht.readlines()
            assert z2 == z1


# DO NOT REMOVE, used to create truth files
def create_files(
    outdir="/tmp/tempsg",
    tg="./tests/python/data/test_genomes",
):
    from pathlib import Path

    from socialgene.cli.nextflow.export_protein_loci_assembly_tables import (
        export_tables,
    )
    from socialgene.config import env_vars

    TEST_GENOMES_DIR = tg
    genome_filepaths = [
        Path(TEST_GENOMES_DIR, i)
        for i in [
            "BGC0000493.gbk",
            "BGC0000492.gbk",
            "lagriamide_mibig_bgc0001946.gbk",
            "GCA_028106715.1_ASM2810671v1_genomic.gbff.gz",
        ]
    ]
    hash_algorithms = ["crc64", "sha512t24u"]
    include_sequences = ["True", "False"]
    for hash_algo in hash_algorithms:
        for include_sequence in include_sequences:
            env_vars["HASHING_ALGORITHM"] = hash_algo
            write_dir = Path(
                Path(outdir, hash_algo, f"seq_included_{include_sequence}")
            )
            write_dir.mkdir(parents=True, exist_ok=True)
            export_tables(
                file_list=genome_filepaths,
                outdir=write_dir,
                n_fasta_splits=1,
                collect_tables_in_memory=True,
                compression=None,
                include_sequences=include_sequence,
            )
