from pathlib import Path
from socialgene.neo4j.sg_modules import (
    SocialgeneModules,
    write_neo4j_headers,
    hmm_sources,
)
import os
import pytest
import itertools
from socialgene.neo4j.sg_modules import hmm_sources as test_hmm_sources
from socialgene.neo4j.admin_import import Neo4jAdminImport
import tempfile
import pathlib


def read_in(fpath):
    if os.path.getsize(fpath) > 1000:
        raise ValueError("File larger than expected")
    with open(fpath, "rb") as f:
        return f.read().decode("utf-8")


expected_headers = {
    "protein_info": "id:ID(protein)\tname\tdescription\tseqlen:int\r\n",
    "amrfinder_hmms_out": ":IGNORE\taccession\tid:ID(amrfinder)\tdescription\tcategory\r\n",
    "local_hmms_out": ":IGNORE\taccession\tid:ID(local)\tdescription\tcategory\r\n",
    "tigrfam_subrole": "id:ID(tigrfam_subrole)\r\n",
    "classiphage_hmms_out": ":IGNORE\taccession\tid:ID(classiphage)\tdescription\tcategory\r\n",
    "tigrfam_mainrole": "id:ID(tigrfam_mainrole)\r\n",
    "protein_to_hmm_header": ":END_ID(protein)\t:START_ID(hmm)\tenv_from:int\tenv_to:int\tseq_pro_score:float\tevalue:float\ti_evalue:float\tdomain_bias:float\tdomain_score:float\tseq_pro_bias:float\thmm_from:int\thmm_to:int\tali_from:int\tali_to:int\r\n",
    "prism_hmms_out": ":IGNORE\taccession\tid:ID(prism)\tdescription\tcategory\r\n",
    "antismash_hmms_out": ":IGNORE\taccession\tid:ID(antismash)\tdescription\tcategory\r\n",
    "amrfinder_hmms_out_relationships": ":START_ID(hmm)\t:IGNORE\t:END_ID(amrfinder)\t:IGNORE\t:IGNORE\r\n",
    "pfam_hmms_out_relationships": ":START_ID(hmm)\t:IGNORE\t:END_ID(pfam)\t:IGNORE\t:IGNORE\r\n",
    "sg_hmm_nodes_out": "id:ID(hmm)\tmodel_length\r\n",
    "mmseqs2": ":START_ID(protein)\t:END_ID(protein)\r\n",
    "blastp": ":START_ID(protein)\t:END_ID(protein)\tpident:float\tlength:int\tmismatch:int\tgapopen:int\tqstart:int\tqend:int\tsstart:int\tsend:int\tevalue:float\tbitscore:float\tqcovhsp:float\r\n",
    "parameters": "id:ID(when)\tSG_LOC_NEO4J\tSG_LOC_HMMS\tNEO4J_dbms_memory_pagecache_size\tNEO4J_dbms_memory_heap_initial__size\tNEO4J_dbms_memory_heap_max__size\tHMMSEARCH_IEVALUE\tHMMSEARCH_BACKGROUND\tHMMSEARCH_BIASFILTER\tHMMSEARCH_NULL2\tHMMSEARCH_SEED\tHMMSEARCH_Z\tHMMSEARCH_DOMZ\tHMMSEARCH_F1\tHMMSEARCH_F2\tHMMSEARCH_F3\tHMMSEARCH_E\tHMMSEARCH_DOME\tHMMSEARCH_INCE\tHMMSEARCH_INCDOME\tHMMSEARCH_BITCUTOFFS\tplatform\tarchitecture\tpy_executable\tpy_version\tgenome_download_command\r\n",
    "locus": "internal_id:ID(nucleotide)\tid\tmol_type\taltitude\tbio_material\tcell_line\tcell_type\tchromosome\tclone\tclone_lib\tcollected_by\tcollection_date\tcountry\tcultivar\tculture_collection\tdb_xref\tdev_stage\tecotype\tenvironmental_sample\tfocus\tgermline\thaplogroup\thaplotype\thost\tidentified_by\tisolate\tisolation_source\tlab_host\tlat_lon\tmacronuclear\tmap\tmating_type\tmetagenome_source\tnote\torganelle\tPCR_primers\tplasmid\tpop_variant\tproviral\trearranged\tsegment\tserotype\tserovar\tsex\tspecimen_voucher\tstrain\tsub_clone\tsubmitter_seqid\tsub_species\tsub_strain\ttissue_lib\ttissue_type\ttransgenic\ttype_material\tvariety\r\n",
    "taxid": "id:ID(taxid)\tname\trank\r\n",
    "mz_cluster_index_nodes": "id:ID(mz_cluster_index)\tcomponent_index:int\tparent_mass:double\tprecursor_mass:double\tsum_precursor_intensity:double\tSmiles:String\trt_mean::double\trt_std_err:double\tlibrary_id:String\tmq_score:double\tmz_error_ppm:double\tmass_diff:String\r\n",
    "resfams_hmms_out": ":IGNORE\taccession\tid:ID(resfams)\tdescription\tcategory\r\n",
    "virus_orthologous_groups_hmms_out": ":IGNORE\taccession\tid:ID(virus_orthologous_groups)\tdescription\tcategory\r\n",
    "tigrfam_hmms_out": ":IGNORE\taccession\tid:ID(tigrfam)\tdescription\tcategory\r\n",
    "assembly_to_taxid": ":START_ID(assembly)\t:END_ID(taxid)\r\n",
    "assembly_to_mz_file": ":START_ID(assembly)\t:END_ID(mz_source_file)\r\n",
    "tigrfamrole_to_mainrole": ":START_ID(tigrfam_role)\t:END_ID(tigrfam_mainrole)\r\n",
    "prism_hmms_out_relationships": ":START_ID(hmm)\t:IGNORE\t:END_ID(prism)\t:IGNORE\t:IGNORE\r\n",
    "assembly_to_locus": ":END_ID(assembly)\t:START_ID(nucleotide)\r\n",
    "resfams_hmms_out_relationships": ":START_ID(hmm)\t:IGNORE\t:END_ID(resfams)\t:IGNORE\t:IGNORE\r\n",
    "classiphage_hmms_out_relationships": ":START_ID(hmm)\t:IGNORE\t:END_ID(classiphage)\t:IGNORE\t:IGNORE\r\n",
    "locus_to_protein": ":START_ID(nucleotide)\t:END_ID(protein)\tstart:int\tend:int\tstrand:int\r\n",
    "bigslice_hmms_out_relationships": ":START_ID(hmm)\t:IGNORE\t:END_ID(bigslice)\t:IGNORE\t:IGNORE\r\n",
    "assembly": "id:ID(assembly)\tmol_type\taltitude\tbio_material\tcell_line\tcell_type\tchromosome\tclone\tclone_lib\tcollected_by\tcollection_date\tcountry\tcultivar\tculture_collection\tdb_xref\tdev_stage\tecotype\tenvironmental_sample\tfocus\tgermline\thaplogroup\thaplotype\thost\tidentified_by\tisolate\tisolation_source\tlab_host\tlat_lon\tmacronuclear\tmap\tmating_type\tmetagenome_source\tnote\torganelle\tPCR_primers\tplasmid\tpop_variant\tproviral\trearranged\tsegment\tserotype\tserovar\tsex\tspecimen_voucher\tstrain\tsub_clone\tsubmitter_seqid\tsub_species\tsub_strain\ttissue_lib\ttissue_type\ttransgenic\ttype_material\tvariety\r\n",
    "tigrfam_to_role": ":START_ID(tigrfam)\t:END_ID(tigrfam_role)\r\n",
    "local_hmms_out_relationships": ":START_ID(hmm)\t:IGNORE\t:END_ID(local)\t:IGNORE\t:IGNORE\r\n",
    "mz_source_file": "id:ID(mz_source_file)\r\n",
    "cluster_to_source_file": ":END_ID(mz_cluster_index)\t:START_ID(mz_source_file)\r\n",
    "bigslice_hmms_out": ":IGNORE\taccession\tid:ID(bigslice)\tdescription\tcategory\r\n",
    "antismash_hmms_out_relationships": ":START_ID(hmm)\t:IGNORE\t:END_ID(antismash)\t:IGNORE\t:IGNORE\r\n",
    "pfam_hmms_out": ":IGNORE\taccession\tid:ID(pfam)\tdescription\tcategory\r\n",
    "tigrfam_to_go": ":START_ID(tigrfam)\t:END_ID(goterm)\r\n",
    "goterm": "id:ID(goterm)\r\n",
    "tigrfam_role": "id:ID(tigrfam_role)\r\n",
    "molecular_network": ":START_ID(mz_cluster_index)\t:END_ID(mz_cluster_index)\tdelta_mz:double\tmeh:float\tcosine:float\tother_score:float\r\n",
    "tigrfamrole_to_subrole": ":START_ID(tigrfam_role)\t:END_ID(tigrfam_subrole)\r\n",
    "tigrfam_hmms_out_relationships": ":START_ID(hmm)\t:IGNORE\t:END_ID(tigrfam)\t:IGNORE\t:IGNORE\r\n",
    "taxid_to_taxid": ":START_ID(taxid)\t:END_ID(taxid)\r\n",
    "virus_orthologous_groups_hmms_out_relationships": ":START_ID(hmm)\t:IGNORE\t:END_ID(virus_orthologous_groups)\t:IGNORE\t:IGNORE\r\n",
}


class Combinator(object):
    instances = list(itertools.combinations(test_hmm_sources, 2)) + list(
        itertools.combinations(test_hmm_sources, 9)
    )


@pytest.fixture(params=Combinator().instances)
def combinator(request):
    return request.param


def test_creation_and_writing_of_neo4j_headers(combinator):
    temp_list = list(SocialgeneModules().nodes.keys())
    temp_list.extend(list(SocialgeneModules().relationships.keys()))
    with tempfile.TemporaryDirectory() as tmpdirname:
        write_neo4j_headers(
            sg_modules=temp_list, hmmlist=hmm_sources, outdir=tmpdirname
        )
        p = Path(tmpdirname).glob("**/*")
        files = [x for x in p if x.is_file()]
        vals = {i.stem: read_in(i) for i in files}
        assert vals == {k: expected_headers[k] for k, v in vals.items()}


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
