from socialgene.neo4j.schema.node_relationship_class import NR


class Relationships(NR):
    def _import(self):
        self.add_relationship(
            neo4j_label="ANNOTATES",
            header_filename="protein_to_hmm_header.header",
            target_subdirectory="parsed_domtblout",
            target_extension="parseddomtblout",
            header=[
                ":END_ID(protein)",
                ":START_ID(hmm)",
                "env_from:Float",
                "env_to:Float",
                "seq_pro_score:Double",
                "evalue:Double",
                "i_evalue:Double",
                "domain_bias:Double",
                "domain_score:Double",
                "seq_pro_bias:Double",
                "hmm_from:Float",
                "hmm_to:Float",
                "ali_from:Float",
                "ali_to:Float",
                "exponentialized:Boolean",
            ],
        )

        self.add_relationship(
            neo4j_label="ASSEMBLES_TO",
            header_filename="assembly_to_locus.header",
            target_subdirectory="genomic_info",
            target_extension="assembly_to_locus",
            header=[":END_ID(assembly)", ":START_ID(nucleotide)"],
        )

        self.add_relationship(
            neo4j_label="ENCODES",
            header_filename="locus_to_protein.header",
            target_subdirectory="genomic_info",
            target_extension="locus_to_protein",
            header=[
                ":START_ID(nucleotide)",
                ":END_ID(protein)",
                "protein_id",
                "locus_tag",
                "start:Float",
                "end:Float",
                "strand:Float",
                "description",
                "partial_on_complete_genome:Boolean",
                "missing_start:Boolean",
                "missing_stop:Boolean",
                "Integerernal_stop:Boolean",
                "partial_in_the_middle_of_a_contig:Boolean",
                "missing_N_terminus:Boolean",
                "missing_C_terminus:Boolean",
                "frameshifted:Boolean",
                "too_short_partial_abutting_assembly_gap:Boolean",
                "incomplete:Boolean",
            ],
        )

        self.add_relationship(
            neo4j_label="TAXON_PARENT",
            header_filename="taxid_to_taxid.header",
            target_subdirectory="taxdump_process",
            target_extension="taxid_to_taxid",
            header=[":START_ID(taxid)", ":END_ID(taxid)"],
        )

        self.add_relationship(
            neo4j_label="GO_ANN",
            header_filename="tigrfam_to_go.header",
            target_subdirectory="tigrfam_info",
            target_extension="tigrfam_to_go",
            header=[":START_ID(hmm_source)", ":END_ID(goterm)"],
        )

        self.add_relationship(
            neo4j_label="PROTEIN_TO_GO",
            header_filename="protein_to_go.header",
            target_subdirectory="protein_info",
            target_extension="protein_to_go",
            header=[":START_ID(protein)", ":END_ID(goterm)"],
        )

        self.add_relationship(
            neo4j_label="GOTERM_RELS",
            multilabel=True,
            header_filename="go_to_go.header",
            target_subdirectory="goterms",
            target_extension="goterm_edgelist",
            header=[":START_ID(goterm)", ":END_ID(goterm)", ":TYPE"],
        )

        self.add_relationship(
            neo4j_label="ROLE_ANN",
            header_filename="tigrfam_to_role.header",
            target_subdirectory="tigrfam_info",
            target_extension="tigrfam_to_role",
            header=[":START_ID(hmm_source)", ":END_ID(tigrfam_role)"],
        )

        self.add_relationship(
            neo4j_label="MAINROLE_ANN",
            header_filename="tigrfamrole_to_mainrole.header",
            target_subdirectory="tigrfam_info",
            target_extension="tigrfamrole_to_mainrole",
            header=[":START_ID(tigrfam_role)", ":END_ID(tigrfam_mainrole)"],
        )

        self.add_relationship(
            neo4j_label="SUBROLE_ANN",
            header_filename="tigrfamrole_to_subrole.header",
            target_subdirectory="tigrfam_info",
            target_extension="tigrfamrole_to_subrole",
            header=[":START_ID(tigrfam_role)", ":END_ID(tigrfam_subrole)"],
        )

        self.add_relationship(
            neo4j_label="IS_TAXON",
            header_filename="assembly_to_taxid.header",
            target_subdirectory="genomic_info",
            target_extension="assembly_to_taxid",
            header=[":START_ID(assembly)", ":END_ID(taxid)"],
        )

        self.add_relationship(
            neo4j_label="BLASTP",
            header_filename="blastp.header",
            target_subdirectory="diamond_blastp",
            target_extension="blast6",
            header=[
                ":START_ID(protein)",
                ":END_ID(protein)",
                "pident:Double",
                "length:Float",
                "mismatch:Float",
                "gapopen:Float",
                "qstart:Float",
                "qend:Float",
                "sstart:Float",
                "send:Float",
                "evalue:Double",
                "bitscore:Double",
                "qcovhsp:Double",
            ],
        )

        self.add_relationship(
            neo4j_label="MMSEQS2",
            header_filename="mmseqs2.header",
            target_subdirectory="mmseqs2_cluster",
            target_extension="mmseqs2_results_cluster.tsv",
            header=[":START_ID(protein)", ":END_ID(protein)"],
        )

        self.add_relationship(
            neo4j_label="CLUSTER_TO_FILE",
            header_filename="cluster_to_source_file.header",
            target_subdirectory="paired_omics",
            target_extension="cluster_to_source_file",
            header=[":END_ID(mz_cluster_index)", ":START_ID(mz_source_file)"],
        )

        self.add_relationship(
            neo4j_label="MOLECULAR_NETWORK",
            header_filename="molecular_network.header",
            target_subdirectory="paired_omics",
            target_extension="molecular_network",
            header=[
                ":START_ID(mz_cluster_index)",
                ":END_ID(mz_cluster_index)",
                "delta_mz:double",
                "meh:Double",
                "cosine:Double",
                "other_score:Double",
            ],
        )

        self.add_relationship(
            neo4j_label="METABO",
            header_filename="assembly_to_mz_file.header",
            target_subdirectory="paired_omics",
            target_extension="assembly_to_mz_file",
            header=[":START_ID(assembly)", ":END_ID(mz_source_file)"],
        )

        self.add_relationship(
            neo4j_label="SOURCE_DB",
            header_filename="hmm_source_relationships.header",
            target_subdirectory="hmm_info",
            target_extension=".hmminfo",
            header=[
                ":END_ID(hmm_source)",
                ":IGNORE",
                ":IGNORE",
                ":IGNORE",
                ":IGNORE",
                ":IGNORE",
                ":IGNORE",
                ":IGNORE",
                ":IGNORE",
                ":START_ID(hmm)",
                ":IGNORE",
                ":IGNORE",
                ":IGNORE",
                ":IGNORE",
                ":IGNORE",
                ":IGNORE",
            ],
        )
