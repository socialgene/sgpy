# import os
# import tempfile
# from pathlib import Path

# import pytest

# from socialgene.cli.export_protein_loci_assembly_tables import export_tables

# FIXTURE_DIR = os.path.dirname(os.path.realpath(__file__))
# FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)
# TEST_GENOMES_DIR = os.path.join(FIXTURE_DIR, "data", "test_genomes")

# EXPECTED_RESULTS_DIR = os.path.join(FIXTURE_DIR, "cli", "expected_results")

# genome_filepaths = [
#     Path(TEST_GENOMES_DIR, i)
#     for i in [
#         "BGC0000493.gbk",
#         "BGC0000492.gbk",
#         "lagriamide_mibig_bgc0001946.gbk",
#         "GCA_028106715.1_ASM2810671v1_genomic.gbff.gz",
#     ]
# ]


# files = [
#     "assembly",
#     "assembly_to_locus",
#     "assembly_to_taxid",
#     "fasta_split_1.faa",
#     "loci",
#     "locus_to_protein",
#     "protein_ids",
#     "protein_info",
# ]

# temp_dir = tempfile.TemporaryDirectory()
# export_tables(
#     file_list=genome_filepaths,
#     outdir=temp_dir.name,
#     n_fasta_splits=1,
#     collect_tables_in_memory=True,
# )


# @pytest.mark.parametrize("item", files)
# def test_gbk_parsing_and_flatfile_for_neo4j_creation_collect_tables_in_memory_false(
#     item,
# ):
#     assert Path(temp_dir.name, item).exists()
#     with open(Path(temp_dir.name, item), "r") as h:
#         z1 = h.readlines()
#         with open(Path(EXPECTED_RESULTS_DIR, item), "r") as ht:
#             z2 = ht.readlines()
#         assert z1 == z2
