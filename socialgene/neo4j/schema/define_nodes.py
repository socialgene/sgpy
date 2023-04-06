from socialgene.base.molbio import SOURCE_KEYS
from socialgene.neo4j.schema.node_relationship_class import NR


class Nodes(NR):
    def _import(self):
        self.add_node(
            neo4j_label="parameters",
            header_filename="parameters.header",
            target_subdirectory="parameters",
            target_extension="socialgene_parameters",
            header=[
                "id:ID(when)",
                "SG_LOC_NEO4J",
                "SG_LOC_HMMS",
                "NEO4J_dbms_memory_pagecache_size",
                "NEO4J_dbms_memory_heap_initial__size",
                "NEO4J_dbms_memory_heap_max__size",
                "HMMSEARCH_IEVALUE",
                "HMMSEARCH_BACKGROUND",
                "HMMSEARCH_BIASFILTER",
                "HMMSEARCH_NULL2",
                "HMMSEARCH_SEED",
                "HMMSEARCH_Z",
                "HMMSEARCH_DOMZ",
                "HMMSEARCH_F1",
                "HMMSEARCH_F2",
                "HMMSEARCH_F3",
                "HMMSEARCH_E",
                "HMMSEARCH_DOME",
                "HMMSEARCH_INCE",
                "HMMSEARCH_INCDOME",
                "HMMSEARCH_BITCUTOFFS",
                "platform",
                "architecture",
                "py_executable",
                "py_version",
                "genome_download_command",
            ],
        )

        self.add_node(
            neo4j_label="assembly",
            header_filename="assembly.header",
            target_subdirectory="genomic_info",
            target_extension="assemblies",
            header=["id:ID(assembly)"] + SOURCE_KEYS,
        )

        self.add_node(
            neo4j_label="nucleotide",
            header_filename="locus.header",
            target_subdirectory="genomic_info",
            target_extension="loci",
            header=["internal_id:ID(nucleotide)"] + ["external_id"] + SOURCE_KEYS,
        )

        self.add_node(
            neo4j_label="protein",
            header_filename="protein_ids.header",
            target_subdirectory="protein_info",
            target_extension="protein_ids",
            header=[
                "id:ID(protein)",
            ],
        )
        self.add_node(
            neo4j_label="protein_access",
            header_filename="protein_info.header",
            target_subdirectory="protein_info",
            target_extension="protein_info",
            header=[
                ":IGNORE",
                "id:ID(protein_access)",
                "description",
                "seqlen:int",
            ],
        )
        self.add_node(
            neo4j_label="hmm",
            header_filename="sg_hmm_nodes.header",
            target_subdirectory="hmm_tsv_parse",
            target_extension="sg_hmm_nodes",
            header=["id:ID(hmm)", "model_length"],
        )

        self.add_node(
            neo4j_label="goterm",
            header_filename="goterm.header",
            target_subdirectory="tigrfam_info",
            target_extension="goterm",
            header=[
                "id:ID(goterm)",
            ],
        )

        self.add_node(
            neo4j_label="tigrfam_role",
            header_filename="tigrfam_role.header",
            target_subdirectory="tigrfam_info",
            target_extension="tigrfam_role",
            header=["id:ID(tigrfam_role)"],
        )
        self.add_node(
            neo4j_label="tigrfam_mainrole",
            header_filename="tigrfam_mainrole.header",
            target_subdirectory="tigrfam_info",
            target_extension="tigrfam_mainrole",
            header=["id:ID(tigrfam_mainrole)"],
        )

        self.add_node(
            neo4j_label="tigrfam_subrole",
            header_filename="tigrfam_subrole.header",
            target_subdirectory="tigrfam_info",
            target_extension="tigrfam_subrole",
            header=["id:ID(tigrfam_subrole)"],
        )

        self.add_node(
            neo4j_label="taxid",
            header_filename="taxid.header",
            target_subdirectory="taxdump_process",
            target_extension="nodes_taxid",
            header=["id:ID(taxid)", "name", "rank"],
        )

        self.add_node(
            neo4j_label="mz_cluster_index",
            header_filename="mz_cluster_index_nodes.header",
            target_subdirectory="paired_omics",
            target_extension="mz_cluster_index_nodes",
            header=[
                "id:ID(mz_cluster_index)",
                "component_index:int",
                "parent_mass:double",
                "precursor_mass:double",
                "sum_precursor_intensity:double",
                "Smiles:String",
                "rt_mean::double",
                "rt_std_err:double",
                "library_id:String",
                "mq_score:double",
                "mz_error_ppm:double",
                "mass_diff:String",
            ],
        )

        self.add_node(
            neo4j_label="mz_source_file",
            header_filename="mz_source_file.header",
            target_subdirectory="paired_omics",
            target_extension="hash.mz_source_file", 
            header=[
                "id:ID(mz_source_file)",
            ],
        )

        self.add_node(
            neo4j_label="hmm_source",
            header_filename="hmm_source.header",
            target_subdirectory="hmm_info",
            target_extension=".hmminfo",
            header= [
                "id:ID(hmm_source)",
                ":LABEL",
                "rel_path:String",
                "name:String",
                "acc:String",
                "notes:String",
                "description:String",
                "date:String",
                "hash:String",
                "hash_used:String",
                "model_length:String,"
                "category:String",
                "subcategory:String",
            ],
        )
                
               