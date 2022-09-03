# import os
# from socialgene.base.socialgene import SocialGene
# import tempfile
# from pathlib import Path
# import pandas as pd
# import pytest
# import hashlib
# from socialgene.entrypoints.export_protein_loci_assembly_tables import export_tables

# FIXTURE_DIR = os.path.dirname(os.path.realpath(__file__))
# FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)
# FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)

# FIXTURE_DIR = os.path.join(FIXTURE_DIR, "data")

# gbk_path = os.path.join(FIXTURE_DIR, "lagriamide_mibig_bgc0001946.gbk")

# sg_object = SocialGene()
# sg_object.parse(gbk_path)
# locus_key = [
#     key
#     for key, value in sg_object.assemblies["lagriamide_mibig_bgc0001946"].loci.items()
#     if "bgc0001946" in key.lower()
# ][0]

# def md6(fname):
#     hash_md5 = hashlib.md5()
#     with open(fname, "rb") as f:
#         for chunk in iter(lambda: f.read(4096), b""):
#             hash_md5.update(chunk)
#     return hash_md5.hexdigest()

# def check_locus_to_protein(tmpdirname):
#     #print(pd.read_csv(Path(tmpdirname, "protein_info"), delimiter="\t", header=None))
#     assert md6(Path(tmpdirname, "locus_to_protein")) == "ae3b20b4f9f1c39657982242ee649438"


# def check_protein_info(tmpdirname):
#     assert hashlib.md5(open(Path(tmpdirname, "protein_info"),'rb').read()).hexdigest() == "a"


# def check_assembly_to_locus(tmpdirname):
#     txt = Path(tmpdirname, "assembly_to_locus").read_text()
#     assert txt == f"lagriamide_mibig_bgc0001946\t{locus_key}\n"


# def check_loci(tmpdirname):
#     txt = Path(tmpdirname, "loci").read_text()
#     temp = locus_key.split("___")[1]
#     assert txt == f"{locus_key}\t{temp}\n"


# def check_assemblies(tmpdirname):
#     txt = Path(tmpdirname, "assemblies").read_text()
#     assert txt == "lagriamide_mibig_bgc0001946\n"


# def test_export_locus_to_protein():

#     with tempfile.TemporaryDirectory() as tmpdirname:
#         sg_object.export_locus_to_protein(outdir=tmpdirname)
#         check_locus_to_protein(tmpdirname=tmpdirname)


# def test_export_protein_info():

#     with tempfile.TemporaryDirectory() as tmpdirname:
#         sg_object.export_protein_info(outdir=tmpdirname)
#         check_protein_info(tmpdirname=tmpdirname)


# def test_export_assembly_to_locus():

#     with tempfile.TemporaryDirectory() as tmpdirname:
#         sg_object.export_assembly_to_locus(outdir=tmpdirname)
#         check_assembly_to_locus(tmpdirname=tmpdirname)


# def test_export_loci():

#     with tempfile.TemporaryDirectory() as tmpdirname:
#         sg_object.export_loci(outdir=tmpdirname)
#         check_loci(tmpdirname=tmpdirname)


# def test_export_assembly():

#     with tempfile.TemporaryDirectory() as tmpdirname:
#         sg_object.export_assembly(outdir=tmpdirname)
#         check_assemblies(tmpdirname=tmpdirname)


# def test_export_tables():

#     glob_path = os.path.join(FIXTURE_DIR, "*mide_mibig_bgc0001946.gbk")

#     with tempfile.TemporaryDirectory() as tmpdirname:
#         export_tables(
#             sequence_files_glob=glob_path, outdir=tmpdirname, n_fasta_splits=5
#         )
#         check_locus_to_protein(tmpdirname=tmpdirname)
#         check_protein_info(tmpdirname=tmpdirname)
#         check_assembly_to_locus(tmpdirname=tmpdirname)
#         check_loci(tmpdirname=tmpdirname)
#         check_assemblies(tmpdirname=tmpdirname)
#         hashes = []
#         for i in range(1, 6):
#             hasher = hashlib.sha256()
#             temp = Path(tmpdirname, f"fasta_split_{i}.faa").read_text()
#             hasher.update(temp.encode("utf-8"))
#             hashes.append(hasher.hexdigest())
#         assert hashes == [
#             "594b75be14dff8f026a6400c0d5ad865edcb1543be6159e571cc35de6bfdcb59",
#             "2708149c3632e6d7599e7b46c546d569040d279fa128b4421c025ba500bdab9e",
#             "7a8f94346dc4cab786e2cd8ffd61d9df5e8ea263755aea5371720a9ceb3b1e1a",
#             "9dd7cf210a6a3aa1ca52d3b1f6869f3727f9d5884c5e11f8cdb60ae49354efc0",
#             "3b8009c0c56a52797c409b6e76903550027fff8b6413e7bbeceff0f9bb9a5e21",
#         ]


# def test_export_tables_fails():
#     glob_path = os.path.join(FIXTURE_DIR, "*mide_mibig_bgc0001946.gbk")

#     with tempfile.TemporaryDirectory() as tmpdirname:
#         with pytest.raises(Exception) as e_info:
#             export_tables(
#                 sequence_files_glob="idontexists", outdir=tmpdirname, n_fasta_splits=5
#             )
#             export_tables(
#                 sequence_files_glob=glob_path, outdir=tmpdirname, n_fasta_splits=5
#             )
#             # neither 0 or 6 should be present
#             with pytest.raises(Exception) as e_info:
#                 Path(tmpdirname, "fasta_split_0.faa").read_text()
#                 print(e_info)  # make linters happy
#             with pytest.raises(Exception) as e_info:
#                 Path(tmpdirname, "fasta_split_6.faa").read_text()
#                 print(e_info)  # make linters happy
