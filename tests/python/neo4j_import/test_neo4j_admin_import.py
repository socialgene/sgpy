from pathlib import Path
import hashlib
from socialgene.neo4j.sg_modules import (
    sgModules,
    write_neo4j_headers,
    hmm_sources,
)


from socialgene.neo4j.admin_import import Neo4jAdminImport
import tempfile
import pathlib


def md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


# @pytest.mark.parametrize(
#     "test_input, expected",
#     [
#         (
#             ["base"],
#             [
#                 "assembly.header",
#                 "assembly_to_locus.header",
#                 "locus.header",
#                 "locus_to_protein.header",
#                 "parameters.header",
#                 "protein_info.header",
#             ],
#         ),
#         (
#             ["base", "hmms", "tigrfam"],
#             [
#                 "amrfinder_hmms_out.header",
#                 "amrfinder_hmms_out_relationships.header",
#                 "antismash_hmms_out.header",
#                 "antismash_hmms_out_relationships.header",
#                 "assembly.header",
#                 "assembly_to_locus.header",
#                 "bigslice_hmms_out.header",
#                 "bigslice_hmms_out_relationships.header",
#                 "classiphage_hmms_out.header",
#                 "classiphage_hmms_out_relationships.header",
#                 "goterm.header",
#                 "local_hmms_out.header",
#                 "local_hmms_out_relationships.header",
#                 "locus.header",
#                 "locus_to_protein.header",
#                 "parameters.header",
#                 "pfam_hmms_out.header",
#                 "pfam_hmms_out_relationships.header",
#                 "prism_hmms_out.header",
#                 "prism_hmms_out_relationships.header",
#                 "protein_info.header",
#                 "protein_to_hmm_header.header",
#                 "resfams_hmms_out.header",
#                 "resfams_hmms_out_relationships.header",
#                 "sg_hmm_nodes_out.header",
#                 "tigrfam_hmms_out.header",
#                 "tigrfam_hmms_out_relationships.header",
#                 "tigrfam_mainrole.header",
#                 "tigrfam_role.header",
#                 "tigrfam_subrole.header",
#                 "tigrfam_to_go.header",
#                 "tigrfam_to_role.header",
#                 "tigrfamrole_to_mainrole.header",
#                 "tigrfamrole_to_subrole.header",
#                 "virus_orthologous_groups_hmms_out.header",
#                 "virus_orthologous_groups_hmms_out_relationships.header",
#             ],
#         ),
#     ],
# )
# def test_filenames(test_input, expected):
#     """Test that parameterized headers are writing"""
#     with tempfile.TemporaryDirectory() as tmpdirname:
#         a = Neo4jImportData(sg_modules=test_input)
#         write_neo4j_headers(sg_modules: list, hmmlist: list, outdir: str)
#         a.write_headers(tmpdirname)
#         p = Path(tmpdirname).glob("**/*")
#         files = [x for x in p if x.is_file()]
#         files.sort()
#         filenames = [i.name for i in files]
#     assert filenames == expected


def test_hashes():
    """Test all the header files' contents"""
    temp_list = list(sgModules().nodes.keys())
    temp_list.extend(list(sgModules().relationships.keys()))
    with tempfile.TemporaryDirectory() as tmpdirname:
        write_neo4j_headers(
            sg_modules=temp_list, hmmlist=hmm_sources, outdir=tmpdirname
        )
        p = Path(tmpdirname).glob("**/*")
        files = [x for x in p if x.is_file()]
        files.sort()
        hashes = {i.name: md5(i) for i in files}
    for k, v in {
        "amrfinder_hmms_out.header": "b031f834168b13799f8a6fb6c82a5638",
        "amrfinder_hmms_out_relationships.header": "aaf4e761a36fba347b54568c8bd29889",
        "antismash_hmms_out.header": "d7277102a630fdc5ccec56d2f4f6ec65",
        "antismash_hmms_out_relationships.header": "3414ec2bbe37b9af520d4a9581da226f",
        "assembly.header": "c9da874fb88e75c70dcd18b3ada3a813",
        "assembly_to_locus.header": "14f9e0e8ffbd25b240dce27e49d9c0ae",
        "assembly_to_mz_file.header": "4400959be4c356ff9ca97ace161293e0",
        "assembly_to_taxid.header": "3f9af852094b60745dc22f1a36c876c8",
        "bigslice_hmms_out.header": "62636ef084bbdc48a3dd7e62f9d0bbe4",
        "bigslice_hmms_out_relationships.header": "30a7cbd8eff47fa11efe7cf2d0bf0e31",
        "blastp.header": "73dcdcb237ea88ed5ccff662110e9b96",
        "classiphage_hmms_out.header": "93e9578f762871457abc220859d4d842",
        "classiphage_hmms_out_relationships.header": "c9662ee87dcf4e041664aaa6b87170a1",
        "cluster_to_source_file.header": "f4de9debcc4dbafc8d5d5b6bf3871d51",
        "goterm.header": "5871ec6e85b2a8947e14f6f1a1b203c7",
        "local_hmms_out.header": "fb7c3134cc56544384329d1db64a3617",
        "local_hmms_out_relationships.header": "11c582c6758a02d7a7b02c5e9f3f4560",
        "locus.header": "6a58a48ab2c8d8f641280aac68d09dfd",
        "locus_to_protein.header": "44e1aae080b0b65a9081d77511efe2c0",
        "mmseqs2.header": "5b191c60314719d8f501c35eee7950c6",
        "molecular_network.header": "4e0095892c523cb462677ff68bb85c69",
        "mz_cluster_index_nodes.header": "4f8f9aa1abcc54a912eb28e94920c15c",
        "mz_source_file.header": "5ac710b20f508bd8c2d87a2c0f1a718f",
        "parameters.header": "bfccc8ca760cdd40b8c6faed661e0d1a",
        "pfam_hmms_out.header": "ef274414c36618474bc9f1d677403e26",
        "pfam_hmms_out_relationships.header": "2fbff8b8d3a5b7fb4c0e7bf2fc8ca7f2",
        "prism_hmms_out.header": "278806347d722935f178973ad8e827f3",
        "prism_hmms_out_relationships.header": "8267a9949d0077ae0a022c042577ad3c",
        "protein_info.header": "bf0197bca389043bd7f5be0bbddab134",
        "protein_to_hmm_header.header": "6a600653c416e901a9f5f57850a7df57",
        "resfams_hmms_out.header": "4ca26515cf1c7eb1360cf74625d0a90c",
        "resfams_hmms_out_relationships.header": "c738437e7779edd3dfc491e06e34158b",
        "sg_hmm_nodes_out.header": "484071a40dcc110a1e7698ea4dbb12b6",
        "taxid.header": "c7e0c717b4cf917544017da85f37c7eb",
        "taxid_to_taxid.header": "cc049bd566ccc90ca64a6cf5ccbd7f5c",
        "tigrfam_hmms_out.header": "502843b111a7d4448a9e16b976ea8fd9",
        "tigrfam_hmms_out_relationships.header": "6baafaf44d2f1b7565fe2d8cd736cb4f",
        "tigrfam_mainrole.header": "128bdbf7acdff546dc046da3da3c3745",
        "tigrfam_role.header": "860575661f5dc766d7bd84709acf4327",
        "tigrfam_subrole.header": "2365fdb76ebf95329c541db265201c63",
        "tigrfam_to_go.header": "fdd5d3e0a4dd4a5e5b2ec1bc3ebdcb27",
        "tigrfam_to_role.header": "bfeb5ce147a7417be63e2adb943474ed",
        "tigrfamrole_to_mainrole.header": "d9c3660b2ace6b35ac9f24db5cd92559",
        "tigrfamrole_to_subrole.header": "597e590650456e37383dd9d061ad2666",
        "virus_orthologous_groups_hmms_out.header": "6ddd3e7662369d3fa548b537cee6f5da",
        "virus_orthologous_groups_hmms_out_relationships.header": "84b210cac16589f2426a51ba27654615",
    }.items():
        assert hashes[k] == v


def test_neo4j_admin_import_dir_creation():
    temp = Neo4jAdminImport(
        neo4j_top_dir="topdir",
        input_sg_modules=["base", "hmms", "mmseqs2"],
        hmmlist=["antismash", "amrfinder"],
        cpus=1,
        additional_args=None,
        uid=100,
        gid=200,
    )
    files = []
    dirs = []
    with tempfile.TemporaryDirectory() as tmpdirname:
        temp.create_neo4j_directories(tmpdirname)
        # iterate directory
        for entry in pathlib.Path(tmpdirname).iterdir():
            # check if it a file
            if entry.is_file():
                files.append(entry.name)
            if entry.is_dir():
                dirs.append(entry.name)
    dirs.sort()
    assert dirs == ["data", "logs", "plugins"]
    assert files == ["import.report"]
