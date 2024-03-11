from socialgene.base.molbio import LocusAssemblyMetadata
from socialgene.neo4j.neo4j_element import Node

# This file can only contain class objects that inherit from Node


class PARAMETERS(Node):
    """Parameters and environmental variables used during database creation"""

    neo4j_label = ["parameters"]
    description = "Parameters and environmental variables used during database creation"
    header_filename = "parameters.header"
    target_subdirectory = "parameters"
    target_extension = "socialgene_parameters"
    header = [
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
    property_specification = {
        "uid": str,
        "SG_LOC_NEO4J": str,
        "SG_LOC_HMMS": str,
        "NEO4J_dbms_memory_pagecache_size": str,
        "NEO4J_dbms_memory_heap_initial__size": str,
        "NEO4J_dbms_memory_heap_max__size": str,
        "HMMSEARCH_IEVALUE": str,
        "HMMSEARCH_BACKGROUND": str,
        "HMMSEARCH_BIASFILTER": str,
        "HMMSEARCH_NULL2": str,
        "HMMSEARCH_SEED": int,
        "HMMSEARCH_Z": int,
        "HMMSEARCH_DOMZ": int,
        "HMMSEARCH_F1": float,
        "HMMSEARCH_F2": float,
        "HMMSEARCH_F3": float,
        "HMMSEARCH_E": float,
        "HMMSEARCH_DOME": float,
        "HMMSEARCH_INCE": float,
        "HMMSEARCH_INCDOME": float,
        "HMMSEARCH_BITCUTOFFS": str,
        "platform": str,
        "architecture": str,
        "py_executable": str,
        "py_version": str,
        "genome_download_command": str,
    }
    constraints_unique = ["uid"]


class ASSEMBLY(Node):
    """Represents a single genome/assembly/BGC. If the input was a FASTA file or if assembly wasn't in the genbank metadata then this will represent the file the data came from."""

    neo4j_label = ["assembly"]
    description = "Represents a single genome/assembly/BGC. If the input was a FASTA file or if assembly wasn't in the genbank metadata then this will represent the file the data came from."
    header_filename = "assembly.header"
    target_subdirectory = "genomic_info"
    target_extension = "assemblies"
    header = ["uid:ID(assembly)"] + sorted(LocusAssemblyMetadata.__slots__)
    constraints_unique = ["uid"]
    property_specification = {"uid": str} | {
        k: str for k in LocusAssemblyMetadata.__slots__
    }
    required_properties = ["uid"]


class NUCLEOTIDE(Node):
    """Represents a single nucleotide sequence (e.g. a contig/scaffold/chromosome)"""

    neo4j_label = ["nucleotide"]
    description = (
        "Represents a single nucleotide sequence (e.g. a contig/scaffold/chromosome)"
    )
    header_filename = "locus.header"
    target_subdirectory = "genomic_info"
    target_extension = "loci"
    header = ["uid:ID(nucleotide)"] + ["external_id"] + LocusAssemblyMetadata.__slots__
    constraints_unique = ["uid"]
    nonunique_index = ["external_id"]
    property_specification = {"uid": str, "external_id": str} | {
        k: str for k in LocusAssemblyMetadata.__slots__
    }
    required_properties = ["uid", "external_id"]


class PROTEIN(Node):
    """Represents a non-redundant protein"""

    neo4j_label = ["protein"]
    description = "Represents a non-redundant protein"
    header_filename = "protein_ids.header"
    target_subdirectory = "protein_info"
    target_extension = "protein_ids"
    constraints_unique = ["uid"]

    def __init__(
        self,
        include_sequences=True,
    ):
        if include_sequences:
            self.header = ["uid:ID(protein)", "crc64", "sequence"]
        else:
            self.header = ["uid:ID(protein)", "crc64"]

    property_specification = {"uid": str, "crc64": str, "sequence": str}


class GOTERM(Node):
    """Represent a GO term"""

    neo4j_label = ["goterm"]
    description = "Represent a GO term"
    header_filename = "goterms.header"
    target_subdirectory = "goterms"
    target_extension = "goterms"
    header = ["uid:ID(goterm)", "name", "namespace", "def"]
    property_specification = {"uid": str, "name": str, "namespace": str}
    constraints_unique = ["uid"]


class TIGRFAM_ROLE(Node):
    """Represents a TIGRFAM role"""

    neo4j_label = ["tigrfam_role"]
    description = "Represents a TIGRFAM role"
    header_filename = "tigrfam_role.header"
    target_subdirectory = "tigrfam_info"
    target_extension = "tigrfam_role"
    header = ["uid:ID(tigrfam_role)"]
    property_specification = {"uid": str}
    constraints_unique = ["uid"]


class TIGRFAM_MAINROLE(Node):
    """Represents a TIGRFAM main role"""

    neo4j_label = ["tigrfam_mainrole"]
    description = "Represents a TIGRFAM main role"
    header_filename = "tigrfam_mainrole.header"
    target_subdirectory = "tigrfam_info"
    target_extension = "tigrfam_mainrole"
    header = ["uid:ID(tigrfam_mainrole)"]
    property_specification = {"uid": str}
    constraints_unique = ["uid"]


class TIGRFAM_SUBROLE(Node):
    """Represents a TIGRFAM sub role"""

    neo4j_label = ["tigrfam_subrole"]
    description = "Represents a TIGRFAM sub role"
    header_filename = "tigrfam_subrole.header"
    target_subdirectory = "tigrfam_info"
    target_extension = "tigrfam_subrole"
    header = ["uid:ID(tigrfam_subrole)"]
    property_specification = {"uid": str}
    constraints_unique = ["uid"]


class TAXID(Node):
    """Represents a single taxon within NCBI taxonomy"""

    neo4j_label = ["taxid"]
    description = "Represents a single taxon within NCBI taxonomy"
    header_filename = "taxid.header"
    target_subdirectory = "taxdump_process"
    target_extension = "nodes_taxid"
    header = ["uid:ID(taxid)", "name", "rank"]
    constraints_unique = ["uid"]
    property_specification = {"uid": str, "name": str, "rank": str}
    required_properties = ["uid"]


class HMM_SOURCE(Node):
    """Represents the source of an HMM model (e.g. PFAM)"""

    neo4j_label = ["hmm_source"]
    description = "Represents the source of an HMM model (e.g. PFAM)"
    header_filename = "hmm_source.header"
    target_subdirectory = "hmm_info"
    target_extension = "hmminfo"
    header = [
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
    property_specification = {
        "uid": str,
        ":LABEL": str,
        "rel_path": str,
        "name": str,
        "acc": str,
        "notes": str,
        "description": str,
        "date": str,
        "hash": str,
        "hash_used": str,
        "model_length": str,
        "super_category": str,
        "category": str,
        "subcategory": str,
        "ga": str,
        "tc": str,
        "nc": str,
    }
    constraints_unique = ["uid"]


class HMM(Node):
    """Represents a single non-redundant HMM model"""

    neo4j_label = ["hmm"]
    description = "Represents a single non-redundant HMM model"
    header_filename = "sg_hmm_nodes.header"
    target_subdirectory = "hmm_info"
    target_extension = "sg_hmm_nodes"
    header = ["uid:ID(hmm)"]
    constraints_unique = ["uid"]
    property_specification = {"uid": str}
    constraints_unique = ["uid"]
