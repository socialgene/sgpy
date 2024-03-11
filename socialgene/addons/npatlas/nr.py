"""https://www.npatlas.org"""

from socialgene.addons.chebi.nr import ChebiNode
from socialgene.addons.chemistry.nr import ChemicalCompoundNode
from socialgene.addons.classyfire.nr import ClassyFireNode
from socialgene.addons.gnps_library.nr import GnpsLibrarySpectrumNode
from socialgene.addons.mibig.nr import Mibig
from socialgene.addons.npclassifier.nr import (
    NPClassifierClass,
    NPClassifierPathway,
    NPClassifierSuperclass,
)
from socialgene.addons.npmrd.nr import Npmrd
from socialgene.addons.publication.nr import Publication
from socialgene.neo4j.neo4j_element import Node, Relationship
from socialgene.nextflow.nodes import TAXID


class NPAtlasNode(Node):
    neo4j_label = ["npatlas"]
    description = "Represents a single NPAtlas entry"
    property_specification = {
        "uid": str,
        "original_name": str,
        "mol_formula": str,
        "mol_weight": float,
        "exact_mass": float,
        "inchikey": str,
        "smiles": str,
        "cluster_id": int,
        "node_id": int,
        "synonyms": str,
        "inchi": str,
        "m_plus_h": float,
        "m_plus_na": float,
        "genus": str,
        "species": str,
    }
    required_properties = ["uid"]
    constraints_unique = ["uid"]


class MibigToNPAtlas(Relationship):
    neo4j_label = "PRODUCES"
    description = "Connects an NPAtlas entry to a Mibig entry"
    start_class = Mibig
    end_class = NPAtlasNode


class NPAtlasToNpmrd(Relationship):
    neo4j_label = "HAS"
    description = "Connects an NPAtlas entry to an Npmrd entry"
    start_class = NPAtlasNode
    end_class = Npmrd


class NPAtlasToGnps(Relationship):
    neo4j_label = "HAS"
    description = "Connects an NPAtlas entry to a GNPS entry"
    start_class = NPAtlasNode
    end_class = GnpsLibrarySpectrumNode


class NPAtlasPublication(Publication):
    pass


class NPAtlasToPublication(Relationship):
    neo4j_label = "HAS"
    description = "Connects an NPAtlas entry to a publication"
    start_class = NPAtlasNode
    end_class = NPAtlasPublication


class NPAtlasToTaxID(Relationship):
    neo4j_label = "PRODUCES"
    description = "Connects an NPAtlas entry to a TaxID entry"
    start_class = TAXID
    end_class = NPAtlasNode


class NPAtlasToClassyFireDirectParent(Relationship):
    neo4j_label = "DIRECT_PARENT"
    description = "Connects an NPAtlas entry to a ClassyFire entry"
    start_class = NPAtlasNode
    end_class = ClassyFireNode


class NPAtlasToClassyFireLowestClass(Relationship):
    neo4j_label = "LOWEST_CLASS"
    description = "Connects an NPAtlas entry to a ClassyFire entry"
    start_class = NPAtlasNode
    end_class = ClassyFireNode


class NPAtlasToClassyFireIntermediateNodes(Relationship):
    neo4j_label = "INTERMEDIATE_NODES"
    description = "Connects an NPAtlas entry to a ClassyFire entry"
    start_class = NPAtlasNode
    end_class = ClassyFireNode


class NPAtlasToClassyFireAlternativeParents(Relationship):
    neo4j_label = "ALTERNATIVE_PARENTS"
    description = "Connects an NPAtlas entry to a ClassyFire entry"
    start_class = NPAtlasNode
    end_class = ClassyFireNode


class NPAtlasToNpclassifierClass(Relationship):
    neo4j_label = "IS_A"
    description = "Connects an NPAtlas entry to a NPClassifierClass entry"
    start_class = NPAtlasNode
    end_class = NPClassifierClass


class NPAtlasToNpclassifierPathway(Relationship):
    neo4j_label = "IS_A"
    description = "Connects an NPAtlas entry to a NPClassifierPathway entry"
    start_class = NPAtlasNode
    end_class = NPClassifierPathway


class NPAtlasToNpclassifierSuperclass(Relationship):
    neo4j_label = "IS_A"
    description = "Connects an NPAtlas entry to a NPClassifierSuperclass entry"
    start_class = NPAtlasNode
    end_class = NPClassifierSuperclass


class NPAtlasToChem(Relationship):
    neo4j_label = "IS_A"
    description = "Connects an NPAtlas entry to a ChemicalCompound entry"
    start_class = NPAtlasNode
    end_class = ChemicalCompoundNode


class NPAtlasToChebi(Relationship):
    neo4j_label = "IS_A"
    description = "Connects an NPAtlas entry to a ChEBI entry"
    start_class = NPAtlasNode
    end_class = ChebiNode
