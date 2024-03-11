import socialgene.nextflow.nodes as Nodes
from socialgene.neo4j.neo4j_element import Relationship

# This file can only contain class objects that inherit from Relationship


class ANNOTATES(Relationship):
    """Represents the relationship between a protein and a HMM"""

    neo4j_label = "ANNOTATES"
    description = "Represents the relationship between a protein and a HMM"
    header_filename = "protein_to_hmm_header.header"
    target_subdirectory = "parsed_domtblout"
    target_extension = "parseddomtblout"
    header = [
        ":END_ID(protein)",
        ":START_ID(hmm)",
        "env_from:Long",
        "env_to:Long",
        "seq_pro_score:Float",
        "evalue:Long",
        "i_evalue:Long",
        "domain_bias:Float",
        "domain_score:Float",
        "seq_pro_bias:Float",
        "hmm_from:Long",
        "hmm_to:Long",
        "ali_from:Long",
        "ali_to:Long",
        "exponentialized:Boolean",
    ]
    start_class = Nodes.HMM
    end_class = Nodes.PROTEIN


class ASSEMBLES_TO(Relationship):
    """Represents the relationship between an assembly and a nucleotide sequence (e.g. scaffold, contig, chromosome)"""

    neo4j_label = "ASSEMBLES_TO"
    description = "Represents the relationship between an assembly and a nucleotide sequence (e.g. scaffold, contig, chromosome)"
    header_filename = "assembly_to_locus.header"
    target_subdirectory = "genomic_info"
    target_extension = "assembly_to_locus"
    header = [":END_ID(assembly)", ":START_ID(nucleotide)"]
    start_class = Nodes.NUCLEOTIDE
    end_class = Nodes.ASSEMBLY


class ENCODES(Relationship):
    """Represents the relationship between a nucleotide sequence and a protein"""

    neo4j_label = "ENCODES"
    description = (
        "Represents the relationship between a nucleotide sequence and a protein"
    )
    header_filename = "locus_to_protein.header"
    target_subdirectory = "genomic_info"
    target_extension = "locus_to_protein"
    header = [
        ":START_ID(nucleotide)",
        ":END_ID(protein)",
        "external_id",
        "locus_tag",
        "start:Long",
        "end:Long",
        "strand:Long",
        "description",
        "partial_on_complete_genome:Boolean",
        "missing_start:Boolean",
        "missing_stop:Boolean",
        "internal_stop:Boolean",
        "partial_in_the_middle_of_a_contig:Boolean",
        "missing_N_terminus:Boolean",
        "missing_C_terminus:Boolean",
        "frameshifted:Boolean",
        "too_short_partial_abutting_assembly_gap:Boolean",
        "incomplete:Boolean",
    ]
    start_class = Nodes.NUCLEOTIDE
    end_class = Nodes.PROTEIN


class TAXON_PARENT(Relationship):
    """Represents the relationship between a taxon and its parent taxon"""

    neo4j_label = "TAXON_PARENT"
    description = "Represents the relationship between a taxon and its parent taxon"
    header_filename = "taxid_to_taxid.header"
    target_subdirectory = "taxdump_process"
    target_extension = "taxid_to_taxid"
    header = [":START_ID(taxid)", ":END_ID(taxid)"]
    start_class = Nodes.TAXID
    end_class = Nodes.TAXID


class GO_ANN(Relationship):
    """Represents the relationship between a TIGRFAM HMM and a GO term"""

    neo4j_label = "GO_ANN"
    description = "Represents the relationship between a TIGRFAM HMM and a GO term"
    header_filename = "tigrfam_to_go.header"
    target_subdirectory = "tigrfam_info"
    target_extension = "tigrfam_to_go"
    header = [":START_ID(hmm_source)", ":END_ID(goterm)"]
    start_class = Nodes.HMM_SOURCE
    end_class = Nodes.GOTERM


class PROTEIN_TO_GO(Relationship):
    """Represents the relationship between a protein and a GO term"""

    neo4j_label = "PROTEIN_TO_GO"
    description = "Represents the relationship between a protein and a GO term"
    header_filename = "protein_to_go.header"
    target_subdirectory = "protein_info"
    target_extension = "protein_to_go"
    header = [":START_ID(protein)", ":END_ID(goterm)"]
    start_class = Nodes.PROTEIN
    end_class = Nodes.GOTERM


class GOTERM_RELS(Relationship):
    """Represents the relationship between GO terms"""

    neo4j_label = "GOTERM_RELS"
    description = "Represents the relationship between GO terms"
    multilabel = True
    header_filename = "go_to_go.header"
    target_subdirectory = "goterms"
    target_extension = "goterm_edgelist"
    header = [":START_ID(goterm)", ":END_ID(goterm)", ":TYPE"]
    start_class = Nodes.GOTERM
    end_class = Nodes.GOTERM


class ROLE_ANN(Relationship):
    """Represents the relationship between a TIGRFAM HMM and a role"""

    neo4j_label = "ROLE_ANN"
    description = "Represents the relationship between a TIGRFAM HMM and a role"
    header_filename = "tigrfam_to_role.header"
    target_subdirectory = "tigrfam_info"
    target_extension = "tigrfam_to_role"
    header = [":START_ID(hmm_source)", ":END_ID(tigrfam_role)"]
    start_class = Nodes.HMM_SOURCE
    end_class = Nodes.TIGRFAM_ROLE


class MAINROLE_ANN(Relationship):
    """Represents the relationship between a TIGRFAM role and a main role"""

    neo4j_label = "MAINROLE_ANN"
    description = "Represents the relationship between a TIGRFAM role and a main role"
    header_filename = "tigrfamrole_to_mainrole.header"
    target_subdirectory = "tigrfam_info"
    target_extension = "tigrfamrole_to_mainrole"
    header = [":START_ID(tigrfam_role)", ":END_ID(tigrfam_mainrole)"]
    start_class = Nodes.TIGRFAM_ROLE
    end_class = Nodes.TIGRFAM_MAINROLE


class SUBROLE_ANN(Relationship):
    """Represents the relationship between a TIGRFAM role and a subrole"""

    neo4j_label = "SUBROLE_ANN"
    description = "Represents the relationship between a TIGRFAM role and a subrole"
    header_filename = "tigrfamrole_to_subrole.header"
    target_subdirectory = "tigrfam_info"
    target_extension = "tigrfamrole_to_subrole"
    header = [":START_ID(tigrfam_role)", ":END_ID(tigrfam_subrole)"]
    start_class = Nodes.TIGRFAM_ROLE
    end_class = Nodes.TIGRFAM_SUBROLE


class IS_TAXON(Relationship):
    """Represents the relationship between an assembly and a taxon"""

    neo4j_label = "IS_TAXON"
    description = "Represents the relationship between an assembly and a taxon"
    header_filename = "assembly_to_taxid.header"
    target_subdirectory = "genomic_info"
    target_extension = "assembly_to_taxid"
    header = [":START_ID(assembly)", ":END_ID(taxid)"]
    start_class = Nodes.ASSEMBLY
    end_class = Nodes.TAXID


class BLASTP(Relationship):
    """Represents the relationship between two proteins calculated with DIAMOND BLASTp"""

    neo4j_label = "BLASTP"
    description = "Represents the relationship between two proteins calculated with DIAMOND BLASTp"
    header_filename = "blastp.header"
    target_subdirectory = "diamond_blastp"
    target_extension = "blast6"
    header = [
        ":START_ID(protein)",
        ":END_ID(protein)",
        "pident:Float",
        "length:Long",
        "mismatch:Long",
        "gapopen:Long",
        "qstart:Long",
        "qend:Long",
        "sstart:Long",
        "send:Long",
        "evalue:Float",
        "bitscore:Float",
        "qcovhsp:Float",
    ]
    start_class = Nodes.PROTEIN
    end_class = Nodes.PROTEIN


class MMSEQS2(Relationship):
    """Represents the relationship between two proteins calculated with MMseqs2"""

    neo4j_label = "MMSEQS2"
    description = (
        "Represents the relationship between two proteins calculated with MMseqs2"
    )
    header_filename = "mmseqs2.header"
    target_subdirectory = "mmseqs2_cluster"
    target_extension = "mmseqs2_results_cluster.tsv"
    header = [":START_ID(protein)", ":END_ID(protein)", ":TYPE"]
    start_class = Nodes.PROTEIN
    end_class = Nodes.PROTEIN


class SOURCE_DB(Relationship):
    """Represents the relationship between a HMM and HMM source database(s)"""

    neo4j_label = "SOURCE_DB"
    description = "Represents the relationship between a HMM and HMM source database(s)"
    header_filename = "hmm_source_relationships.header"
    target_subdirectory = "hmm_info"
    target_extension = ".hmminfo"
    header = [
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
        ":IGNORE",
    ]
    start_class = Nodes.HMM
    end_class = Nodes.HMM_SOURCE
