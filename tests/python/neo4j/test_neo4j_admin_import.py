import os
import pathlib
import tempfile
from pathlib import Path

import pytest

from socialgene.neo4j.schema.socialgene_modules import SocialgeneModules
from socialgene.neo4j.utils.admin_import import Neo4jAdminImport


def read_in(fpath):
    if os.path.getsize(fpath) > 1000:
        raise ValueError("File larger than expected")
    with open(fpath, "rb") as f:
        return f.read().decode("utf-8")


expected_headers = {
    "tigrfam_subrole": "uid:ID(tigrfam_subrole)\r\n",
    "tigrfam_mainrole": "uid:ID(tigrfam_mainrole)\r\n",
    "protein_to_hmm_header": ":END_ID(protein)\t:START_ID(hmm)\tenv_from:Long\tenv_to:Long\tseq_pro_score:Float\tevalue:Long\ti_evalue:Long\tdomain_bias:Float\tdomain_score:Float\tseq_pro_bias:Float\thmm_from:Long\thmm_to:Long\tali_from:Long\tali_to:Long\texponentialized:Boolean\r\n",
    "mmseqs2": ":START_ID(protein)\t:END_ID(protein)\t:TYPE\r\n",
    "blastp": ":START_ID(protein)\t:END_ID(protein)\tpident:Float\tlength:Long\tmismatch:Long\tgapopen:Long\tqstart:Long\tqend:Long\tsstart:Long\tsend:Long\tevalue:Float\tbitscore:Float\tqcovhsp:Float\r\n",
    "parameters": "uid:ID(when)\tSG_LOC_NEO4J\tSG_LOC_HMMS\tNEO4J_dbms_memory_pagecache_size\tNEO4J_dbms_memory_heap_initial__size\tNEO4J_dbms_memory_heap_max__size\tHMMSEARCH_IEVALUE\tHMMSEARCH_BACKGROUND\tHMMSEARCH_BIASFILTER\tHMMSEARCH_NULL2\tHMMSEARCH_SEED\tHMMSEARCH_Z\tHMMSEARCH_DOMZ\tHMMSEARCH_F1\tHMMSEARCH_F2\tHMMSEARCH_F3\tHMMSEARCH_E\tHMMSEARCH_DOME\tHMMSEARCH_INCE\tHMMSEARCH_INCDOME\tHMMSEARCH_BITCUTOFFS\tplatform\tarchitecture\tpy_executable\tpy_version\tgenome_download_command\r\n",
    "locus": "uid:ID(nucleotide)\texternal_id\taltitude\tbio_material\tbioproject\tbiosample\tcell_line\tcell_type\tchromosome\tclone\tclone_lib\tcollected_by\tcollection_date\tcountry\tcultivar\tculture_collection\tdb_xref\tdev_stage\tecotype\tenvironmental_sample\tfocus\tgermline\thaplogroup\thaplotype\thost\tidentified_by\tisolate\tisolation_source\tlab_host\tlat_lon\tmacronuclear\tmap\tmating_type\tmetagenome_source\tmol_type\tnote\torganelle\torganism\tpcr_primers\tplasmid\tpop_variant\tproviral\trearranged\tsegment\tserotype\tserovar\tsex\tspecimen_voucher\tstrain\tsub_clone\tsubmitter_seqid\tsub_species\tsub_strain\ttissue_lib\ttissue_type\ttransgenic\ttype_material\tvariety\r\n",
    "taxid": "uid:ID(taxid)\tname\trank\r\n",
    "sg_hmm_nodes": "uid:ID(hmm)\r\n",
    "assembly_to_taxid": ":START_ID(assembly)\t:END_ID(taxid)\r\n",
    "tigrfamrole_to_mainrole": ":START_ID(tigrfam_role)\t:END_ID(tigrfam_mainrole)\r\n",
    "assembly_to_locus": ":END_ID(assembly)\t:START_ID(nucleotide)\r\n",
    "locus_to_protein": ":START_ID(nucleotide)\t:END_ID(protein)\texternal_id\tlocus_tag\tstart:Long\tend:Long\tstrand:Long\tdescription\tpartial_on_complete_genome:Boolean\tmissing_start:Boolean\tmissing_stop:Boolean\tinternal_stop:Boolean\tpartial_in_the_middle_of_a_contig:Boolean\tmissing_N_terminus:Boolean\tmissing_C_terminus:Boolean\tframeshifted:Boolean\ttoo_short_partial_abutting_assembly_gap:Boolean\tincomplete:Boolean\r\n",
    "assembly": "uid:ID(assembly)\taltitude\tbio_material\tbioproject\tbiosample\tcell_line\tcell_type\tchromosome\tclone\tclone_lib\tcollected_by\tcollection_date\tcountry\tcultivar\tculture_collection\tdb_xref\tdev_stage\tecotype\tenvironmental_sample\tfocus\tgermline\thaplogroup\thaplotype\thost\tidentified_by\tisolate\tisolation_source\tlab_host\tlat_lon\tmacronuclear\tmap\tmating_type\tmetagenome_source\tmol_type\tnote\torganelle\torganism\tpcr_primers\tplasmid\tpop_variant\tproviral\trearranged\tsegment\tserotype\tserovar\tsex\tspecimen_voucher\tstrain\tsub_clone\tsub_species\tsub_strain\tsubmitter_seqid\ttissue_lib\ttissue_type\ttransgenic\ttype_material\tvariety\r\n",
    "tigrfam_to_role": ":START_ID(hmm_source)\t:END_ID(tigrfam_role)\r\n",
    "protein_to_go": ":START_ID(protein)\t:END_ID(goterm)\r\n",
    "protein_ids": "uid:ID(protein)\tcrc64\tsequence\r\n",
    "hmm_source_relationships": ":END_ID(hmm_source)\t:IGNORE\t:IGNORE\t:IGNORE\t:IGNORE\t:IGNORE\t:IGNORE\t:IGNORE\t:IGNORE\t:START_ID(hmm)\t:IGNORE\t:IGNORE\t:IGNORE\t:IGNORE\t:IGNORE\t:IGNORE\t:IGNORE\r\n",
    "go_to_go": ":START_ID(goterm)\t:END_ID(goterm)\t:TYPE\r\n",
    "tigrfam_to_go": ":START_ID(hmm_source)\t:END_ID(goterm)\r\n",
    "tigrfam_role": "uid:ID(tigrfam_role)\r\n",
    "hmm_source": "uid:ID(hmm_source)\t:LABEL\trel_path:String\tname:String\tacc:String\tnotes:String\tdescription:String\tdate:String\thash:String\thash_used:String\tmodel_length:String\tsuper_category:String\tcategory:String\tsubcategory:String\tga:String\ttc:String\tnc:String\r\n",
    "tigrfamrole_to_subrole": ":START_ID(tigrfam_role)\t:END_ID(tigrfam_subrole)\r\n",
    "taxid_to_taxid": ":START_ID(taxid)\t:END_ID(taxid)\r\n",
    "goterms": "uid:ID(goterm)\tname\tnamespace\tdef\r\n",
}


@pytest.fixture(scope="session")
def tempor_dir(tmpdir_factory):
    sg_mod = SocialgeneModules()
    sg_mod.add_modules(list(sg_mod.modules.keys()))
    fn = tmpdir_factory.mktemp("data")
    sg_mod.write_neo4j_headers(outdir=fn)
    return fn


@pytest.mark.parametrize("k", [k for k in expected_headers.keys()])
def test_creation_and_writing_of_neo4j_headers(tempor_dir, k):
    p = Path(tempor_dir).glob("**/*")
    files = [x for x in p if x.is_file()]
    vals = {i.stem: read_in(i) for i in files}
    assert vals[k] == expected_headers[k]


def test_neo4j_admin_import_dir_creation():
    temp = Neo4jAdminImport(
        neo4j_top_dir="topdir",
        module_list=["base", "base_hmm", "mmseqs"],
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
