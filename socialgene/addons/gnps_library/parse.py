import re
from typing import List

from rdkit import Chem

from socialgene.addons.base import ExternalBaseClass
from socialgene.addons.npclassifier.nr import NPClassifierClass
from socialgene.base.chem import ChemicalCompound
from socialgene.neo4j.neo4j_element import Node, Relationship
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
