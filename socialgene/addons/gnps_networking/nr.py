from uuid import uuid4
from socialgene.addons.gnps_library.nr import GnpsLibrarySpectrumNode


from socialgene.neo4j.neo4j_element import Node, Relationship
from socialgene.nextflow.nodes import ASSEMBLY
from socialgene.utils.logging import log


class ClusterNode(Node):
    neo4j_label = "gnps_cluster"
    description = "Represents a GNPS molecular networking cluster"
    required_properties = ["cluster_index", "workflow_uuid"]
    property_specification = {
        "workflow_uuid": str,
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
    }





class SpectrumNode(Node):
    neo4j_label = "spectrum"
    description = "Represents a GNPS molecular networking spectrum"
    property_specification = {
        "uid": str,
        "original_filename": str,
        "parentmass": float,
        "charge": int,
        "rettime": float,
        "assembly": str,
    }


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
    required_properties=[]


class ClusterToAssembly(Relationship):
    neo4j_label = "FROM"
    description = "Connects a GNPS cluster to a GNPS assembly"
    start_class = ClusterNode
    end_class = ASSEMBLY
    property_specification = {
    }
    required_properties = []
