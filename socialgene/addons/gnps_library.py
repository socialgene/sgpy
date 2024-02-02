import re
from typing import List

from socialgene.addons.base import ExternalBaseClass

from socialgene.addons.npclassifier import NPClassifierClass
from socialgene.neo4j.neo4j import GraphDriver
from socialgene.neo4j.neo4j_element import Node, Relationship
from socialgene.utils.logging import log


class GnpsLibrarySpectrumNode(Node):
    def __init__(self, *args, **kwargs):
        super().__init__(
            neo4j_label="gnps_library_spectrum",
            description="Represents a GNPS library spectrum",
            properties={
                "uid": "string",
                "compound_name": "string",
                "compound_source": "string",
                "pi": "string",
                "data_collector": "string",
                "adduct": "string",
                "precursor_mz": "float",
                "exactmass": "float",
                "charge": "int",
                "cas_number": "string",
                "pubmed_id": "string",
                "smiles": "string",
                "inchi": "string",
                "inchi_aux": "string",
                "library_class": "string",
                "ionmode": "string",
                "libraryqualitystring": "string",
                "mqscore": "float",
                "tic_query": "float",
                "rt_query": "float",
                "mzerrorppm": "float",
                "sharedpeaks": "int",
                "massdiff": "float",
                "libmz": "float",
                "specmz": "float",
                "speccharge": "int",
                "moleculeexplorerdatasets": "string",
                "moleculeexplorerfiles": "string",
                "molecular_formula": "string",
                "inchikey": "string",
                "inchikey_planar": "string",
            },
        )


class IonSourceNode(Node):
    def __init__(self, *args, **kwargs):
        super().__init__(
            neo4j_label="ion_source",
            description="Represents an ion source",
            properties={
                "uid": "string",
            },
        )


class InstrumentNode(Node):
    def __init__(self, *args, **kwargs):
        super().__init__(
            neo4j_label="instrument",
            description="Represents an instrument",
            properties={
                "uid": "string",
            },
        )


class OrganismNode(Node):
    def __init__(self, *args, **kwargs):
        super().__init__(
            neo4j_label="organism",
            description="Represents an organism (as defined by GNPS)",
            properties={
                "uid": "string",
            },
        )


class FromIonRel(Relationship):
    def __init__(self, *args, **kwargs):
        super().__init__(
            neo4j_label="FROM",
            description="Connects a GNPS spectrum to an ion source",
            start=GnpsLibrarySpectrumNode,
            end=IonSourceNode,
        )


class FromInstrumentRel(Relationship):
    def __init__(self, *args, **kwargs):
        super().__init__(
            neo4j_label="FROM",
            description="Connects a GNPS spectrum to an instrument source",
            start=GnpsLibrarySpectrumNode,
            end=InstrumentNode,
        )


class FromOrganismRel(Relationship):
    def __init__(self, *args, **kwargs):
        super().__init__(
            neo4j_label="FROM",
            description="Connects a GNPS spectrum to an organism (as defined by GNPS)",
            start=GnpsLibrarySpectrumNode,
            end=OrganismNode,
        )


class GnpsLibrarySpectrumIsA(Relationship):
    def __init__(self, *args, **kwargs):
        super().__init__(
            neo4j_label="IS_A",
            description="Represents a relationship between gnps_library_spectrum nodes",
            start=GnpsLibrarySpectrumNode,
            end=NPClassifierClass,
        )


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
