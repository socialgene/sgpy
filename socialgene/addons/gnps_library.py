import re
from typing import List

from rdkit import Chem

from socialgene.addons.base import ExternalBaseClass
from socialgene.addons.npclassifier import NPClassifierClass
from socialgene.base.chem import ChemicalCompound
from socialgene.neo4j.neo4j_element import Node, Relationship
from socialgene.utils.logging import log


class GnpsLibrarySpectrumNode(Node):
    neo4j_label = "gnps_library_spectrum"
    description = "Represents a GNPS library spectrum"
    uid = ["uid"]
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


class IonSourceNode(Node):
    neo4j_label = "ion_source"
    description = "Represents an ion source"
    property_specification = {
        "uid": str,
    }


class InstrumentNode(Node):
    neo4j_label = "instrument"
    description = "Represents an instrument"
    property_specification = {
        "uid": str,
    }


class OrganismNode(Node):
    neo4j_label = "organism"
    description = "Represents an organism (as defined by GNPS)"
    property_specification = {
        "uid": str,
    }


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


class GnpsLibrarySpectrumIsA(Relationship):
    neo4j_label = "IS_A"
    description = "Represents a relationship between gnps_library_spectrum nodes"
    start_class = GnpsLibrarySpectrumNode
    end_class = NPClassifierClass


class GnpsLibrarySpectrum(ExternalBaseClass):
    __slots__ = [
        "uid",
        "compound_name",
        "compound_source",
        "pi",
        "data_collector",
        "adduct",
        "precursor_mz",
        "exactmass",
        "charge",
        "cas_number",
        "pubmed_id",
        "smiles",
        "inchi",
        "inchi_aux",
        "library_class",
        "ionmode",
        "libraryqualitystring",
        "mqscore",
        "tic_query",
        "rt_query",
        "mzerrorppm",
        "sharedpeaks",
        "massdiff",
        "libmz",
        "specmz",
        "speccharge",
        "moleculeexplorerdatasets",
        "moleculeexplorerfiles",
        "molecular_formula",
        "inchikey",
        "inchikey_planar",
    ]

    def add_cmpd(self):
        temp = None
        try:
            inchi = self.inchi
            inchi = inchi.removeprefix('"').removesuffix('"')
            temp = ChemicalCompound(self.inchi)
        except Exception:
            pass
        if not temp:
            try:
                temp = ChemicalCompound(self.smiles)
            except Exception:
                pass
        if temp:
            temp.add_to_neo4j()
            self._add_to_neo4j(
                """
                    with $uid as row and $temp as temp
                    MATCH (gls:gnps_library_spectrum {uid: row.uid})
                    MATCH (cc:chemical_compound {uid: temp.uid})
                    MERGE (gls)-[:IS_A]->(cc)
                        """,
                uid=self.uid,
                temp=(Chem.MolToInchi(temp.mol), Chem.MolToSmiles(temp.mol)),
            )

    def __init__(self, **kwargs) -> None:
        for i in self.__slots__:
            setattr(self, i, None)
        for k, v in kwargs.items():
            setattr(self, k, v)

    @staticmethod
    def _extract_CCMSLIB(x) -> List:
        return re.findall("CCMSLIB[0-9]{11}", x)

    @property
    def not_null_properties(self):
        return {
            i: getattr(self, i) for i in self.__slots__ if getattr(self, i) is not None
        }

    def add_node_to_neo4j(self):
        self._add_to_neo4j(
            """
            with $row as row
            MERGE (gls:gnps_library_spectrum {uid: row.uid})
                     SET gls += row
                        """,
            row=self.not_null_properties,
        )

    def link_to_chemical_node(self):
        cmpd = None
        try:
            cmpd = ChemicalCompound(self.inchikey)
        except Exception:
            pass
        if not cmpd:
            try:
                cmpd = ChemicalCompound(self.smiles)
            except Exception:
                pass
        if cmpd:
            log.warning(f"Could not parse a chemical compound for {self.uid}")
            return

        self._add_to_neo4j(
            """
            with $row as row
            MERGE (gls:gnps_library_spectrum {uid: row.uid})
                     SET gls += row
                        """,
            row=self.not_null_properties,
        )
