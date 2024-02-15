"""https://www.npatlas.org"""

import json
from pathlib import Path
import tempfile
from itertools import batched

from socialgene.addons.base import ExternalBaseClass
from socialgene.addons.gnps_library.gnps_library import GnpsLibrarySpectrumNode
from socialgene.addons.mibig.mibig import Mibig
from socialgene.addons.npmrd.npmrd import Npmrd
from socialgene.addons.publication.publication import Publication
from socialgene.base.chem import ChemicalCompound
from socialgene.neo4j.neo4j_element import Node, Relationship
from socialgene.utils.download import download as downloader
from socialgene.utils.logging import log

NPATALAS_URL = "https://www.npatlas.org/static/downloads/NPAtlas_download.json"


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


class NPAtlasParser:
    __slots__ = [
        "entry",
        "uid",
        "original_name",
        "mol_formula",
        "mol_weight",
        "exact_mass",
        "inchikey",
        "smiles",
        "cluster_id",
        "node_id",
        "synonyms",
        "inchi",
        "m_plus_h",
        "m_plus_na",
        "origin_reference",
        "npclassifier",
        "external_ids",
        "ncbi_taxid",
        "genus",
        "species",
        "gnps",
        "npmrd",
        "mibig",
        "classyfire_class",
        "classyfire_subclass",
    ]

    def __init__(self, entry) -> None:
        self.entry = entry
        self._parse_single_entry()

    def _parse_single_entry(self):
        self.ncbi_taxid = None
        self.genus = None
        self.species = None
        self.gnps = set()
        self.npmrd = set()
        self.mibig = set()
        self.classyfire_class = None
        self.classyfire_subclass = None
        self.uid = self.entry.get("npaid", None)
        self.original_name = self.entry.get("original_name", None)
        self.mol_formula = self.entry.get("mol_formula", None)
        self.mol_weight = float(self.entry.get("mol_weight", None))
        self.exact_mass = float(self.entry.get("exact_mass", None))
        self.inchikey = str(self.entry.get("inchikey", None))
        self.smiles = str(self.entry.get("smiles", None))
        self.cluster_id = int(self.entry.get("cluster_id", None))
        self.node_id = int(self.entry.get("node_id", None))
        self.inchi = self.entry.get("inchi", None)
        self.m_plus_h = self.entry.get("m_plus_h", None)
        self.m_plus_na = self.entry.get("m_plus_na", None)
        ######
        if self.entry.get("origin_reference", None):
            self.origin_reference = NPAtlasPublication(
                **self.entry.get("origin_reference")
            )
        # self.synonyms = ";".join(self.entry.get("synonyms", None))
        self.synonyms = None
        self._assign_to_taxon()
        self._assign_to_classy()
        self._assign_external_ids()
        # self.external_ids = self.entry.get("external_ids", None)

    def _assign_to_classy(self):
        try:
            self.classyfire_class = self.entry["classyfire"]["class"][
                "chemont_id"
            ].removeprefix("CHEMONTID:")
        except Exception as e:
            log.debug(e)
        try:
            self.classyfire_subclass = self.entry["classyfire"]["subclass"][
                "chemont_id"
            ].removeprefix("CHEMONTID:")
        except Exception as e:
            log.debug(e)

    def _assign_to_taxon(self):
        try:
            self.ncbi_taxid = self.entry["origin_organism"]["taxon"]["ncbi_id"]
        except Exception:
            pass
        try:
            self.genus = self.entry["origin_organism"]["genus"]
        except Exception:
            pass
        try:
            self.species = self.entry["origin_organism"]["species"]
        except Exception:
            pass

    def _assign_external_ids(self):
        try:
            for external_id in self.entry["external_ids"]:
                match external_id["external_db_name"]:
                    case "gnps":
                        gnps_ids = GnpsLibrarySpectrumNode._extract_CCMSLIB(
                            external_id["external_db_code"]
                        )
                        for uid in gnps_ids:
                            self.gnps.add(GnpsLibrarySpectrumNode(properties={"uid":uid}))
                    case "npmrd":
                        self.npmrd.add(Npmrd(properties={"uid":external_id["external_db_code"]}))
                    case "mibig":
                        self.mibig.add(Mibig(properties={"uid":external_id["external_db_code"]}))
                    case _:
                        pass
        except Exception as e:
            log.debug(e)




class NPAtlas(ExternalBaseClass):
    def __init__(self, url=NPATALAS_URL, atlas_json_path=None) -> None:
        super().__init__()
        self.entries = []
        self.path = atlas_json_path
        if not Path(atlas_json_path).exists():
            self.path=self._download_npatlas(url=url, outpath=atlas_json_path)

    def _download_npatlas(self, outpath, url=NPATALAS_URL):
        if Path(outpath).exists():
            log.debug(f"File already exists at {outpath}")
        else:
            downloader(url, outpath)





    def merge_with_classyfire(self):
        if self.classyfire_class:
            self._add_to_neo4j(
                """
                MERGE (np:npatlas {uid: $uid})
                MERGE (bgc:chemont {uid: $chemont})
                MERGE (bgc)<-[:IS_CLASS]-(np)
                """,
                uid=self.uid,
                chemont=self.classyfire_class,
            )
        if self.classyfire_subclass:
            self._add_to_neo4j(
                """
                MERGE (np:npatlas {uid: $uid})
                MERGE (bgc:chemont {uid: $chemont})
                MERGE (bgc)<-[:IS_SUBCLASS]-(np)
                """,
                uid=self.uid,
                chemont=self.classyfire_subclass,
            )

