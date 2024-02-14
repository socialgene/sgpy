from socialgene.base.molbio import LocusAssemblyMetadata
from socialgene.neo4j.neo4j_element import Node

# This file can only contain class objects that inherit from Node


class PARAMETERS(Node):
    """Parameters and environmental variables used during database creation"""

    neo4j_label="parameters"
    description="Parameters and environmental variables used during database creation"
    header_filename="parameters.header"
    target_subdirectory="parameters"
    target_extension="socialgene_parameters"
    header=[
        "uid:ID(when)",
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
    ]



class ASSEMBLY(Node):
    """Represents a single genome/assembly/BGC. If the input was a FASTA file or if assembly wasn't in the genbank metadata then this will represent the file the data came from."""


    neo4j_label="assembly"
    description="Represents a single genome/assembly/BGC. If the input was a FASTA file or if assembly wasn't in the genbank metadata then this will represent the file the data came from."
    header_filename="assembly.header"
    target_subdirectory="genomic_info"
    target_extension="assemblies"
    header=["uid:ID(assembly)"] + sorted(LocusAssemblyMetadata.__slots__)
    unique_constraints=["uid"]


class NUCLEOTIDE(Node):
    """Represents a single nucleotide sequence (e.g. a contig/scaffold/chromosome)"""


    neo4j_label="nucleotide"
    description="Represents a single nucleotide sequence (e.g. a contig/scaffold/chromosome)"
    header_filename="locus.header"
    target_subdirectory="genomic_info"
    target_extension="loci"
    header=["uid:ID(nucleotide)"] + ["external_id"] + LocusAssemblyMetadata.__slots__
    unique_constraints=["uid"]
    nonunique_index=["external_id"]



class PROTEIN(Node):
    """Represents a non-redundant protein"""
    neo4j_label="protein"
    description="Represents a non-redundant protein"
    header_filename="protein_ids.header"
    target_subdirectory="protein_info"
    target_extension="protein_ids"
    unique_constraints=["uid"]

    def __init__(self, include_sequences=True,):
        if include_sequences:
            self.header = ["uid:ID(protein)", "crc64", "sequence"]
        else:
            self.header = ["uid:ID(protein)", "crc64"]




class GOTERM(Node):
    """Represent a GO term"""
    neo4j_label="goterm"
    description="Represent a GO term"
    header_filename="goterms.header"
    target_subdirectory="goterms"
    target_extension="goterms"
    header=["uid:ID(goterm)", "name", "namespace", "def"]



class TIGRFAM_ROLE(Node):
    """Represents a TIGRFAM role"""


    neo4j_label="tigrfam_role"
    description="Represents a TIGRFAM role"
    header_filename="tigrfam_role.header"
    target_subdirectory="tigrfam_info"
    target_extension="tigrfam_role"
    header=["uid:ID(tigrfam_role)"]



class TIGRFAM_MAINROLE(Node):
    """Represents a TIGRFAM main role"""


    neo4j_label="tigrfam_mainrole"
    description="Represents a TIGRFAM main role"
    header_filename="tigrfam_mainrole.header"
    target_subdirectory="tigrfam_info"
    target_extension="tigrfam_mainrole"
    header=["uid:ID(tigrfam_mainrole)"]



class TIGRFAM_SUBROLE(Node):
    """Represents a TIGRFAM sub role"""

    neo4j_label="tigrfam_subrole"
    description="Represents a TIGRFAM sub role"
    header_filename="tigrfam_subrole.header"
    target_subdirectory="tigrfam_info"
    target_extension="tigrfam_subrole"
    header=["uid:ID(tigrfam_subrole)"]



class TAXID(Node):
    """Represents a single taxon within NCBI taxonomy"""


    neo4j_label="taxid"
    description="Represents a single taxon within NCBI taxonomy"
    header_filename="taxid.header"
    target_subdirectory="taxdump_process"
    target_extension="nodes_taxid"
    header=["uid:ID(taxid)", "name", "rank"]
    unique_constraints=["uid"]



class HMM_SOURCE(Node):
    """Represents the source of an HMM model (e.g. PFAM)"""


    neo4j_label="hmm_source"
    description="Represents the source of an HMM model (e.g. PFAM)"
    header_filename="hmm_source.header"
    target_subdirectory="hmm_info"
    target_extension="hmminfo"
    header=[
        "uid:ID(hmm_source)",
        ":LABEL",
        "rel_path:String",
        "name:String",
        "acc:String",
        "notes:String",
        "description:String",
        "date:String",
        "hash:String",
        "hash_used:String",
        "model_length:String",
        "super_category:String",
        "category:String",
        "subcategory:String",
        "ga:String",
        "tc:String",
        "nc:String",
    ]



class HMM(Node):
    """Represents a single non-redundant HMM model"""


    neo4j_label="hmm"
    description="Represents a single non-redundant HMM model"
    header_filename="sg_hmm_nodes.header"
    target_subdirectory="hmm_info"
    target_extension="sg_hmm_nodes"
    header=["uid:ID(hmm)"]
    unique_constraints=["uid"]

