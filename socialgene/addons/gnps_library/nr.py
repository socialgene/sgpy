import re
from typing import List

from socialgene.addons.chemistry.nr import ChemicalCompoundNode
from socialgene.addons.npclassifier.nr import (
    NPClassifierClass,
    NPClassifierPathway,
    NPClassifierSuperclass,
)
from socialgene.neo4j.neo4j_element import Node, Relationship


class GnpsLibrarySpectrumNode(Node):
    neo4j_label = ["gnps_library_spectrum"]
    description = "Represents a GNPS library spectrum"
    uid = ["uid"]
    constraints_unique = ["uid"]
    required_properties = ["uid"]
    property_specification = {
        "uid": str,
        "compound_name": str,
        "compound_source": str,
        "pi": str,
        "data_collector": str,
        "adduct": str,
        "precursor_mz": float,
        "exactmass": float,
        "charge": int,
        "cas_number": str,
        "pubmed_id": str,
        "smiles": str,
        "inchi": str,
        "inchi_aux": str,
        "library_class": str,
        "ionmode": str,
        "libraryqualitystring": str,
        "mqscore": float,
        "tic_query": float,
        "rt_query": float,
        "mzerrorppm": float,
        "sharedpeaks": int,
        "massdiff": float,
        "libmz": float,
        "specmz": float,
        "speccharge": int,
        "moleculeexplorerdatasets": str,
        "moleculeexplorerfiles": str,
        "molecular_formula": str,
        "inchikey": str,
        "inchikey_planar": str,
    }

    @staticmethod
    def _extract_all_CCMSLIB(x) -> List:
        return re.findall("CCMSLIB[0-9]{11}", x)


class IonSourceNode(Node):
    neo4j_label = ["ion_source"]
    description = "Represents an ion source"
    property_specification = {
        "uid": str,
    }
    constraints_unique = ["uid"]


class InstrumentNode(Node):
    neo4j_label = ["instrument"]
    description = "Represents an instrument"
    property_specification = {
        "uid": str,
    }
    constraints_unique = ["uid"]


class OrganismNode(Node):
    neo4j_label = ["gnps_organism"]
    description = "Represents an organism (as defined by GNPS)"
    property_specification = {
        "uid": str,
    }
    constraints_unique = ["uid"]


class FromIonRel(Relationship):
    neo4j_label = "FROM"
    description = "Connects a GNPS spectrum to an ion source"
    start_class = GnpsLibrarySpectrumNode
    end_class = IonSourceNode


class FromInstrumentRel(Relationship):
    neo4j_label = "FROM"
    description = "Connects a GNPS spectrum to an instrument source"
    start_class = GnpsLibrarySpectrumNode
    end_class = InstrumentNode


class FromOrganismRel(Relationship):
    neo4j_label = "FROM"
    description = "Connects a GNPS spectrum to an organism (as defined by GNPS)"
    start_class = GnpsLibrarySpectrumNode
    end_class = OrganismNode


class GnpsLibrarySpectrumToNPClassifierClass(Relationship):
    neo4j_label = "IS_A"
    description = "Represents a relationship between gnps_library_spectrum nodes"
    start_class = GnpsLibrarySpectrumNode
    end_class = NPClassifierClass


class GnpsLibrarySpectrumToNPClassifierSuperclass(Relationship):
    neo4j_label = "IS_A"
    description = "Represents a relationship between gnps_library_spectrum nodes"
    start_class = GnpsLibrarySpectrumNode
    end_class = NPClassifierSuperclass


class GnpsLibrarySpectrumToNPClassifierPathway(Relationship):
    neo4j_label = "IS_A"
    description = "Represents a relationship between gnps_library_spectrum nodes"
    start_class = GnpsLibrarySpectrumNode
    end_class = NPClassifierPathway


class GnpsLibraryToChem(Relationship):
    neo4j_label = "IS_A"
    description = "Connects a GNPS library spectrum to a chemical compound"
    start_class = GnpsLibrarySpectrumNode
    end_class = ChemicalCompoundNode
