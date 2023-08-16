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
    "protein_info",
]
hash_algo = ["crc64", "sha512t24u"]  # , "sha256"]


@pytest.mark.parametrize("item", files)
@pytest.mark.parametrize("hash_algo", hash_algo)
def test_gbk_parsing_and_flatfile_for_neo4j_creation_collect_tables_in_memory_false(
    item, hash_algo
):
    env_vars["HASHING_ALGORITHM"] = hash_algo
    with tempfile.TemporaryDirectory() as tmpdirname:
        export_tables(
            file_list=genome_filepaths,
            outdir=tmpdirname,
            n_fasta_splits=1,
            collect_tables_in_memory=True,
            compression=None,
        )
        assert Path(tmpdirname, item).exists()
        with open(Path(tmpdirname, item), "r") as h:
            z1 = h.readlines()
            with open(Path(TRUTHS_DIR, hash_algo, item), "r") as ht:
                z2 = ht.readlines()
            assert z1 == z2


def create_files(
    outdir="/tmp/tempsg",
    tg="/home/chase/Documents/github/kwan_lab/socialgene/sgpy/tests/python/data/test_genomes",
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
    hash_algo = ["crc64", "sha512t24u"]  # , "sha256"]
    for i in hash_algo:
        env_vars["HASHING_ALGORITHM"] = i
        Path(outdir, i).mkdir(parents=True, exist_ok=True)
        export_tables(
            file_list=genome_filepaths,
            outdir=Path(outdir, i),
            n_fasta_splits=1,
            collect_tables_in_memory=True,
            compression=None,
        )
