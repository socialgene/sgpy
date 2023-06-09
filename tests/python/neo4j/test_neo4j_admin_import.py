import os
from pathlib import Path
import tempfile
import pathlib

from socialgene.neo4j.schema.socialgene_modules import SocialgeneModules
from socialgene.neo4j.admin_import import Neo4jAdminImport


def read_in(fpath):
    if os.path.getsize(fpath) > 1000:
        raise ValueError("File larger than expected")
    with open(fpath, "rb") as f:
        return f.read().decode("utf-8")


expected_headers = {
    "protein_to_hmm_header": ":END_ID(protein)\t:START_ID(hmm)\tenv_from:int\tenv_to:int\tseq_pro_score:float\tevalue:float\ti_evalue:float\tdomain_bias:float\tdomain_score:float\tseq_pro_bias:float\thmm_from:int\thmm_to:int\tali_from:int\tali_to:int\r\n",
    "mmseqs2": ":START_ID(protein)\t:END_ID(protein)\r\n",
    "blastp": ":START_ID(protein)\t:END_ID(protein)\tpident:float\tlength:int\tmismatch:int\tgapopen:int\tqstart:int\tqend:int\tsstart:int\tsend:int\tevalue:float\tbitscore:float\tqcovhsp:float\r\n",
    "parameters": "uid:ID(when)\tSG_LOC_NEO4J\tSG_LOC_HMMS\tNEO4J_dbms_memory_pagecache_size\tNEO4J_dbms_memory_heap_initial__size\tNEO4J_dbms_memory_heap_max__size\tHMMSEARCH_IEVALUE\tHMMSEARCH_BACKGROUND\tHMMSEARCH_BIASFILTER\tHMMSEARCH_NULL2\tHMMSEARCH_SEED\tHMMSEARCH_Z\tHMMSEARCH_DOMZ\tHMMSEARCH_F1\tHMMSEARCH_F2\tHMMSEARCH_F3\tHMMSEARCH_E\tHMMSEARCH_DOME\tHMMSEARCH_INCE\tHMMSEARCH_INCDOME\tHMMSEARCH_BITCUTOFFS\tplatform\tarchitecture\tpy_executable\tpy_version\tgenome_download_command\r\n",
    "locus": "internal_id:ID(nucleotide)\texternal_id\tmol_type\taltitude\tbio_material\tbioproject\tbiosample\tcell_line\tcell_type\tchromosome\tclone\tclone_lib\tcollected_by\tcollection_date\tcountry\tcultivar\tculture_collection\tdb_xref\tdev_stage\tecotype\tenvironmental_sample\tfocus\tgermline\thaplogroup\thaplotype\thost\tidentified_by\tisolate\tisolation_source\tlab_host\tlat_lon\tmacronuclear\tmap\tmating_type\tmetagenome_source\tnote\torganelle\tPCR_primers\tplasmid\tpop_variant\tproviral\trearranged\tsegment\tserotype\tserovar\tsex\tspecimen_voucher\tstrain\tsub_clone\tsubmitter_seqid\tsub_species\tsub_strain\ttissue_lib\ttissue_type\ttransgenic\ttype_material\tvariety\r\n",
    "taxid": "uid:ID(taxid)\tname\trank\r\n",
    "mz_cluster_index_nodes": "uid:ID(mz_cluster_index)\tcomponent_index:int\tparent_mass:double\tprecursor_mass:double\tsum_precursor_intensity:double\tSmiles:String\trt_mean::double\trt_std_err:double\tlibrary_id:String\tmq_score:double\tmz_error_ppm:double\tmass_diff:String\r\n",
    "sg_hmm_nodes": "uid:ID(hmm)\r\n",
    "assembly_to_taxid": ":START_ID(assembly)\t:END_ID(taxid)\r\n",
    "assembly_to_mz_file": ":START_ID(assembly)\t:END_ID(mz_source_file)\r\n",
    "assembly_to_locus": ":END_ID(assembly)\t:START_ID(nucleotide)\r\n",
    "locus_to_protein": ":START_ID(nucleotide)\t:END_ID(protein)\tlocus_tag\tstart:int\tend:int\tstrand:int\tdescription\tpartial_on_complete_genome:boolean\tmissing_start:boolean\tmissing_stop:boolean\tinternal_stop:boolean\tpartial_in_the_middle_of_a_contig:boolean\tmissing_N_terminus:boolean\tmissing_C_terminus:boolean\tframeshifted:boolean\ttoo_short_partial_abutting_assembly_gap:boolean\tincomplete:boolean\r\n",
    "assembly": "uid:ID(assembly)\tmol_type\taltitude\tbio_material\tbioproject\tbiosample\tcell_line\tcell_type\tchromosome\tclone\tclone_lib\tcollected_by\tcollection_date\tcountry\tcultivar\tculture_collection\tdb_xref\tdev_stage\tecotype\tenvironmental_sample\tfocus\tgermline\thaplogroup\thaplotype\thost\tidentified_by\tisolate\tisolation_source\tlab_host\tlat_lon\tmacronuclear\tmap\tmating_type\tmetagenome_source\tnote\torganelle\tPCR_primers\tplasmid\tpop_variant\tproviral\trearranged\tsegment\tserotype\tserovar\tsex\tspecimen_voucher\tstrain\tsub_clone\tsubmitter_seqid\tsub_species\tsub_strain\ttissue_lib\ttissue_type\ttransgenic\ttype_material\tvariety\r\n",
    "protein_ids": "uid:ID(protein)\r\n",
    "mz_source_file": "uid:ID(mz_source_file)\r\n",
    "hmm_source_relationships": ":END_ID(hmm_source)\t:IGNORE\t:IGNORE\t:IGNORE\t:IGNORE\t:IGNORE\t:IGNORE\t:IGNORE\t:IGNORE\t:START_ID(hmm)\t:IGNORE\t:IGNORE\t:IGNORE\t:IGNORE\t:IGNORE\t:IGNORE\r\n",
    "cluster_to_source_file": ":END_ID(mz_cluster_index)\t:START_ID(mz_source_file)\r\n",
    "hmm_source": "uid:ID(hmm_source)\t:LABEL\trel_path:String\tname:String\tacc:String\tnotes:String\tdescription:String\tdate:String\thash:String\thash_used:String\tmodel_length:String\tcategory:String\tsubcategory:String\tga:String\ttc:String\tnc:String\r\n",
    "molecular_network": ":START_ID(mz_cluster_index)\t:END_ID(mz_cluster_index)\tdelta_mz:double\tmeh:float\tcosine:float\tother_score:float\r\n",
    "taxid_to_taxid": ":START_ID(taxid)\t:END_ID(taxid)\r\n",
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
