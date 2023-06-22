import os
import pathlib
import tempfile
from pathlib import Path

from socialgene.neo4j.admin_import import Neo4jAdminImport
from socialgene.neo4j.schema.socialgene_modules import SocialgeneModules


def read_in(fpath):
    if os.path.getsize(fpath) > 1000:
        raise ValueError("File larger than expected")
    with open(fpath, "rb") as f:
        return f.read().decode("utf-8")


expected_headers = {
    "tigrfam_subrole": "uid:ID(tigrfam_subrole)\r\n",
    "tigrfam_mainrole": "uid:ID(tigrfam_mainrole)\r\n",
    "protein_to_hmm_header": ":END_ID(protein)\t:START_ID(hmm)\tenv_from:Long\tenv_to:Long\tseq_pro_score:Float\tevalue:Long\ti_evalue:Long\tdomain_bias:Float\tdomain_score:Float\tseq_pro_bias:Float\thmm_from:Long\thmm_to:Long\tali_from:Long\tali_to:Long\texponentialized:Boolean\r\n",
    "mmseqs2": ":START_ID(protein)\t:END_ID(protein)\r\n",
    "blastp": ":START_ID(protein)\t:END_ID(protein)\tpident:Double\tlength:Long\tmismatch:Long\tgapopen:Long\tqstart:Long\tqend:Long\tsstart:Long\tsend:Long\tevalue:Double\tbitscore:Double\tqcovhsp:Double\r\n",
    "parameters": "uid:ID(when)\tSG_LOC_NEO4J\tSG_LOC_HMMS\tNEO4J_dbms_memory_pagecache_size\tNEO4J_dbms_memory_heap_initial__size\tNEO4J_dbms_memory_heap_max__size\tHMMSEARCH_IEVALUE\tHMMSEARCH_BACKGROUND\tHMMSEARCH_BIASFILTER\tHMMSEARCH_NULL2\tHMMSEARCH_SEED\tHMMSEARCH_Z\tHMMSEARCH_DOMZ\tHMMSEARCH_F1\tHMMSEARCH_F2\tHMMSEARCH_F3\tHMMSEARCH_E\tHMMSEARCH_DOME\tHMMSEARCH_INCE\tHMMSEARCH_INCDOME\tHMMSEARCH_BITCUTOFFS\tplatform\tarchitecture\tpy_executable\tpy_version\tgenome_download_command\r\n",
    "locus": "uid:ID(nucleotide)\texternal_id\tmol_type\taltitude\tbio_material\tbioproject\tbiosample\tcell_line\tcell_type\tchromosome\tclone\tclone_lib\tcollected_by\tcollection_date\tcountry\tcultivar\tculture_collection\tdb_xref\tdev_stage\tecotype\tenvironmental_sample\tfocus\tgermline\thaplogroup\thaplotype\thost\tidentified_by\tisolate\tisolation_source\tlab_host\tlat_lon\tmacronuclear\tmap\tmating_type\tmetagenome_source\tnote\torganelle\tPCR_primers\tplasmid\tpop_variant\tproviral\trearranged\tsegment\tserotype\tserovar\tsex\tspecimen_voucher\tstrain\tsub_clone\tsubmitter_seqid\tsub_species\tsub_strain\ttissue_lib\ttissue_type\ttransgenic\ttype_material\tvariety\r\n",
    "taxid": "uid:ID(taxid)\tname\trank\r\n",
    "mz_cluster_index_nodes": "uid:ID(mz_cluster_index)\tcomponent_index:Long\tparent_mass:Float\tprecursor_mass:Float\tsum_precursor_intensity:Float\tSmiles:String\trt_mean::Float\trt_std_err:Float\tlibrary_id:String\tmq_score:Float\tmz_error_ppm:Float\tmass_diff:String\r\n",
    "sg_hmm_nodes": "uid:ID(hmm)\r\n",
    "assembly_to_taxid": ":START_ID(assembly)\t:END_ID(taxid)\r\n",
    "assembly_to_mz_file": ":START_ID(assembly)\t:END_ID(mz_source_file)\r\n",
    "tigrfamrole_to_mainrole": ":START_ID(tigrfam_role)\t:END_ID(tigrfam_mainrole)\r\n",
    "assembly_to_locus": ":END_ID(assembly)\t:START_ID(nucleotide)\r\n",
    "locus_to_protein": ":START_ID(nucleotide)\t:END_ID(protein)\tprotein_id\tlocus_tag\tstart:Long\tend:Long\tstrand:Long\tdescription\tpartial_on_complete_genome:Boolean\tmissing_start:Boolean\tmissing_stop:Boolean\tinternal_stop:Boolean\tpartial_in_the_middle_of_a_contig:Boolean\tmissing_N_terminus:Boolean\tmissing_C_terminus:Boolean\tframeshifted:Boolean\ttoo_short_partial_abutting_assembly_gap:Boolean\tincomplete:Boolean\r\n",
    "assembly": "uid:ID(assembly)\tmol_type\taltitude\tbio_material\tbioproject\tbiosample\tcell_line\tcell_type\tchromosome\tclone\tclone_lib\tcollected_by\tcollection_date\tcountry\tcultivar\tculture_collection\tdb_xref\tdev_stage\tecotype\tenvironmental_sample\tfocus\tgermline\thaplogroup\thaplotype\thost\tidentified_by\tisolate\tisolation_source\tlab_host\tlat_lon\tmacronuclear\tmap\tmating_type\tmetagenome_source\tnote\torganelle\tPCR_primers\tplasmid\tpop_variant\tproviral\trearranged\tsegment\tserotype\tserovar\tsex\tspecimen_voucher\tstrain\tsub_clone\tsubmitter_seqid\tsub_species\tsub_strain\ttissue_lib\ttissue_type\ttransgenic\ttype_material\tvariety\r\n",
    "tigrfam_to_role": ":START_ID(hmm_source)\t:END_ID(tigrfam_role)\r\n",
    "protein_to_go": ":START_ID(protein)\t:END_ID(goterm)\r\n",
    "protein_ids": "uid:ID(protein)\r\n",
    "mz_source_file": "uid:ID(mz_source_file)\r\n",
    "hmm_source_relationships": ":END_ID(hmm_source)\t:IGNORE\t:IGNORE\t:IGNORE\t:IGNORE\t:IGNORE\t:IGNORE\t:IGNORE\t:IGNORE\t:START_ID(hmm)\t:IGNORE\t:IGNORE\t:IGNORE\t:IGNORE\t:IGNORE\t:IGNORE\r\n",
    "cluster_to_source_file": ":END_ID(mz_cluster_index)\t:START_ID(mz_source_file)\r\n",
    "go_to_go": ":START_ID(goterm)\t:END_ID(goterm)\t:TYPE\r\n",
    "tigrfam_to_go": ":START_ID(hmm_source)\t:END_ID(goterm)\r\n",
    "tigrfam_role": "uid:ID(tigrfam_role)\r\n",
    "hmm_source": "uid:ID(hmm_source)\t:LABEL\trel_path:String\tname:String\tacc:String\tnotes:String\tdescription:String\tdate:String\thash:String\thash_used:String\tmodel_length:String\tcategory:String\tsubcategory:String\tga:String\ttc:String\tnc:String\r\n",
    "molecular_network": ":START_ID(mz_cluster_index)\t:END_ID(mz_cluster_index)\tdelta_mz:double\tmeh:Double\tcosine:Double\tother_score:Double\r\n",
    "tigrfamrole_to_subrole": ":START_ID(tigrfam_role)\t:END_ID(tigrfam_subrole)\r\n",
    "taxid_to_taxid": ":START_ID(taxid)\t:END_ID(taxid)\r\n",
    "goterms": "uid:ID(goterm)\tname\tnamespace\tdef\r\n",
}


def test_creation_and_writing_of_neo4j_headers():
    sg_mod = SocialgeneModules()
    sg_mod.add_modules(sg_mod.modules.keys())
    with tempfile.TemporaryDirectory() as tmpdirname:
        sg_mod.write_neo4j_headers(outdir=tmpdirname)
        p = Path(tmpdirname).glob("**/*")
        files = [x for x in p if x.is_file()]
        vals = {i.stem: read_in(i) for i in files}
        assert vals == {k: expected_headers[k] for k, v in vals.items()}


def test_neo4j_admin_import_dir_creation():
    temp = Neo4jAdminImport(
        neo4j_top_dir="topdir",
        module_list=["base", "hmms", "mmseqs2"],
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
