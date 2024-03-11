from socialgene.addons.gnps_library.nr import GnpsLibrarySpectrumNode
from socialgene.neo4j.neo4j_element import Node, Relationship
from socialgene.nextflow.nodes import ASSEMBLY


class ClusterNode(Node):
    neo4j_label = ["gnps_cluster"]
    description = "Represents a GNPS molecular networking cluster"
    required_properties = ["cluster_index", "workflow_uuid", "task"]
    property_specification = {
        "task": str,
        "defaultgroups": str,
        "g1": int,
        "g2": int,
        "g3": int,
        "g4": int,
        "g5": int,
        "g6": int,
        "gnpslinkout_cluster": str,
        "gnpslinkout_network": str,
        "mqscore": float,
        "mzerrorppm": float,
        "massdiff": float,
        "rtmean": float,
        "rtmean_min": float,
        "rtstderr": float,
        "uniquefilesources": str,
        "uniquefilesourcescount": int,
        "cluster_index": int,
        "componentindex": int,
        "number_of_spectra": int,
        "parent_mass": float,
        "precursor_charge": int,
        "precursor_mass": float,
        "sumprecursor_intensity": float,
        "workflow_uuid": str,
    }
    constraints_unique = ["cluster_index", "workflow_uuid", "task"]


class SpectrumNode(Node):
    neo4j_label = ["ms2_spectrum"]
    description = "Represents a GNPS molecular networking spectrum"
    property_specification = {
        "specidx": str,
        "original_filename": str,
        "parentmass": float,
        "charge": int,
        "rettime": float,
        "workflow_uuid": str,
    }
    constraints_unique = ["original_filename", "specidx", "workflow_uuid"]


class MassSpecFileNode(Node):
    neo4j_label = ["mass_spectrum_file"]
    description = "Represents a GNPS molecular networking spectrum file"
    property_specification = {
        "filename": str,
        "gnps_filename": str,
        "workflow_uuid": str,
    }
    constraints_unique = ["original_filename", "workflow_uuid"]


class MassSpecFileToSpectrum(Relationship):
    neo4j_label = "HAS"
    description = "Connects a GNPS mass spectrum file to a GNPS spectrum"
    start_class = MassSpecFileNode
    end_class = SpectrumNode
    property_specification = {}
    required_properties = []


class MassSpecFileToAssembly(Relationship):
    neo4j_label = "ANALYSIS_OF"
    description = "Connects a GNPS mass spectrum file to an assembly"
    start_class = MassSpecFileNode
    end_class = ASSEMBLY
    property_specification = {}
    required_properties = []


class ClusterToSpectrum(Relationship):
    neo4j_label = "CLUSTERS_TO"
    description = "Connects a GNPS cluster to a GNPS spectrum"
    start_class = SpectrumNode
    end_class = ClusterNode
    property_specification = {}
    required_properties = []


class LibraryHitRel(Relationship):
    neo4j_label = "LIBRARY_HIT"
    description = "Connects a GNPS cluster to a GNPS library hit"
    start_class = ClusterNode
    end_class = GnpsLibrarySpectrumNode


class MolecularNetwork(Relationship):
    neo4j_label = "MOLECULAR_NETWORK"
    description = "Connects the GNPS molecular networks"
    start_class = ClusterNode
    end_class = ClusterNode
    property_specification = {
        "delta_mz": float,
        "meh": float,
        "cosine": float,
        "otherscore": float,
        "componentindex": int,
        "edgeannotation": str,
    }
    required_properties = []
