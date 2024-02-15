"""https://www.npatlas.org"""

from socialgene.addons.gnps_library.nr import GnpsLibrarySpectrumNode
from socialgene.addons.mibig.nr import Mibig
from socialgene.addons.npmrd.nr import Npmrd
from socialgene.addons.publication.nr import Publication
from socialgene.base.chem import ChemicalCompound
from socialgene.neo4j.neo4j_element import Node, Relationship
from socialgene.utils.download import download as downloader


class NPAtlasNode(Node):
        neo4j_label = "npatlas"
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
            # "origin_reference": Publication,
            "npclassifier": str,
            "external_ids": str,
            # "ncbi_taxid": str,
            "genus": str,
            "species": str,
            # "gnps": GnpsLibrarySpectrum,
            # "npmrd": Npmrd,
            # "mibig": Mibig,
            "classyfire_class": str,
            "classyfire_subclass": str,
        }
        required_properties = ["uid"]


class NPAtlasToMibig(Relationship):
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

class NPAtlasToChemont(Relationship):
    neo4j_label = "IS_CLASS"
    description = "Connects an NPAtlas entry to a Chemont entry"
    start_class = NPAtlasNode
    end_class = ChemicalCompound

class NPAtlasToSubclass(Relationship):
    neo4j_label = "IS_SUBCLASS"
    description = "Connects an NPAtlas entry to a Chemont entry"
    start_class = NPAtlasNode
    end_class = ChemicalCompound

