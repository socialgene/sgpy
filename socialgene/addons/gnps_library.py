import re
from typing import List

from socialgene.addons.base import ExternalBaseClass
from socialgene.neo4j.neo4j import GraphDriver
from socialgene.utils.logging import log


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
        for k, v in kwargs.items():
            setattr(self, k, v)

    @staticmethod
    def _extract_CCMSLIB(x) -> List:
        return re.findall("CCMSLIB[0-9]{11}", x)

    @property
    def props(self):
        return {i: getattr(self, i) for i in self.__slots__}
