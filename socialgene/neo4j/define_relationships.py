from socialgene.neo4j.module_abc import Neo4j_Module

# See Modules.add_node() for a description of what's going on here
# Add a node by copying and modifying a Modules().add_node() below


class Relationships:
    relationships = list()

    def _add(self, **kwargs):
        self.relationships.append(Neo4j_Module(**kwargs))


Relationships()._add(
    neo4j_label="ANNOTATES",
    header_filename="protein_to_hmm_header.header",
    target_subdirectory="parsed_domtblout",
    target_extension="parseddomtblout",
    header=[
        ":END_ID(protein)",
        ":START_ID(hmm)",
        "env_from:int",
        "env_to:int",
        "seq_pro_score:float",
        "evalue:float",
        "i_evalue:float",
        "domain_bias:float",
        "domain_score:float",
        "seq_pro_bias:float",
        "hmm_from:int",
        "hmm_to:int",
        "ali_from:int",
        "ali_to:int",
    ],
)

Relationships()._add(
    neo4j_label="ASSEMBLES_TO",
    header_filename="assembly_to_locus.header",
    target_subdirectory="genomic_info",
    target_extension="assembly_to_locus",
    header=[":END_ID(assembly)", ":START_ID(nucleotide)"],
)

Relationships()._add(
    neo4j_label="ENCODES",
    header_filename="locus_to_protein.header",
    target_subdirectory="genomic_info",
    target_extension="locus_to_protein",
    header=[
        ":START_ID(nucleotide)",
        ":END_ID(protein)",
        "start:int",
        "end:int",
        "strand:int",
    ],
)

Relationships()._add(
    neo4j_label="BELONGS_TO",
    header_filename="taxid_to_taxid.header",
    target_subdirectory="taxdump_process",
    target_extension="taxid_to_taxid",
    header=[":START_ID(taxid)", ":END_ID(taxid)"],
)

Relationships()._add(
    neo4j_label="GO_ANN",
    header_filename="tigrfam_to_go.header",
    target_subdirectory="tigrfam_info",
    target_extension="tigrfam_to_go",
    header=[":START_ID(tigrfam)", ":END_ID(goterm)"],
)

Relationships()._add(
    neo4j_label="ROLE_ANN",
    header_filename="tigrfam_to_role.header",
    target_subdirectory="tigrfam_info",
    target_extension="tigrfam_to_role",
    header=[":START_ID(tigrfam)", ":END_ID(tigrfam_role)"],
)

Relationships()._add(
    neo4j_label="MAINROLE_ANN",
    header_filename="tigrfamrole_to_mainrole.header",
    target_subdirectory="tigrfam_info",
    target_extension="tigrfamrole_to_mainrole",
    header=[":START_ID(tigrfam_role)", ":END_ID(tigrfam_mainrole)"],
)

Relationships()._add(
    neo4j_label="SUBROLE_ANN",
    header_filename="tigrfamrole_to_subrole.header",
    target_subdirectory="tigrfam_info",
    target_extension="tigrfamrole_to_subrole",
    header=[":START_ID(tigrfam_role)", ":END_ID(tigrfam_subrole)"],
)

Relationships()._add(
    neo4j_label="TAXONOMY",
    header_filename="assembly_to_taxid.header",
    target_subdirectory="genomic_info",
    target_extension="assembly_to_taxid",
    header=[":START_ID(assembly)", ":END_ID(taxid)"],
)

Relationships()._add(
    neo4j_label="BLASTP",
    header_filename="blastp.header",
    target_subdirectory="diamond_blastp",
    target_extension="blast6",
    header=[
        ":START_ID(protein)",
        ":END_ID(protein)",
        "pident:float",
        "length:int",
        "mismatch:int",
        "gapopen:int",
        "qstart:int",
        "qend:int",
        "sstart:int",
        "send:int",
        "evalue:float",
        "bitscore:float",
        "qcovhsp:float",
    ],
)

Relationships()._add(
    neo4j_label="MMSEQS2",
    header_filename="mmseqs2.header",
    target_subdirectory="mmseqs2_easycluster",
    target_extension="mmseqs2_results_cluster.tsv",
    header=[":START_ID(protein)", ":END_ID(protein)"],
)

Relationships()._add(
    neo4j_label="CLUSTER_TO_FILE",
    header_filename="cluster_to_source_file.header",
    target_subdirectory="paired_omics",
    target_extension="cluster_to_source_file",
    header=[":END_ID(mz_cluster_index)", ":START_ID(mz_source_file)"],
)

Relationships()._add(
    neo4j_label="MOLECULAR_NETWORK",
    header_filename="molecular_network.header",
    target_subdirectory="paired_omics",
    target_extension="molecular_network",
    header=[
        ":START_ID(mz_cluster_index)",
        ":END_ID(mz_cluster_index)",
        "delta_mz:double",
        "meh:float",
        "cosine:float",
        "other_score:float",
    ],
)

Relationships()._add(
    neo4j_label="METABO",
    header_filename="assembly_to_mz_file.header",
    target_subdirectory="paired_omics",
    target_extension="assembly_to_mz_file",
    header=[":START_ID(assembly)", ":END_ID(mz_source_file)"],
)
