"""https://www.npatlas.org"""
import json
import tempfile

import requests
from rich.progress import Progress

from socialgene.addons.base import ExternalBaseClass
from socialgene.addons.gnps_library import GnpsLibrarySpectrum
from socialgene.addons.mibig import Mibig
from socialgene.addons.npmrd import Npmrd
from socialgene.addons.publication import Publication
from socialgene.base.chem import ChemicalCompound
from socialgene.utils.download import download as downloader
from socialgene.utils.logging import log
from itertools import batched
NPATALAS_URL = "https://www.npatlas.org/static/downloads/NPAtlas_download.json"

def download(url=NPATALAS_URL, outpath=None):
    if outpath:
        downloader(url, outpath)
    else:
        with tempfile.NamedTemporaryFile() as tf:
            downloader(url, tf.name)

class NPAtlasPublication(Publication):
    def __init__(self, doi, pmid, authors, title, journal, year, **kwargs) -> None:
        super().__init__()
        self.doi = self._extract_doi(doi)
        self.pmid = str(pmid)
        self.authors = str(authors)
        self.title = str(title)
        self.journal = str(journal)
        self.year = str(year)

# class NPAtlasNode(Node):
#     ...
class NPAtlasEntry:
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
            self.origin_reference = NPAtlasPublication(**self.entry.get("origin_reference"))
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
                        gnps_ids = GnpsLibrarySpectrum._extract_CCMSLIB(
                            external_id["external_db_code"]
                        )
                        for uid in gnps_ids:
                            self.gnps.add(GnpsLibrarySpectrum(uid=uid))
                    case "npmrd":
                        self.npmrd.add(Npmrd(uid=external_id["external_db_code"]))
                    case "mibig":
                        self.mibig.add(Mibig(uid=external_id["external_db_code"]))
                    case _:
                        pass
        except Exception as e:
            log.debug(e)
    @property
    def _node_prop_dict(self):
        return {
            "uid": self.uid,
            "original_name": self.original_name,
            "mol_formula": self.mol_formula,
            "mol_weight": self.mol_weight,
            "exact_mass": self.exact_mass,
            "inchikey": self.inchikey,
            "smiles": self.smiles,
            "cluster_id": self.cluster_id,
            "node_id": self.node_id,
            "synonyms": self.synonyms,
            "inchi": self.inchi,
            "m_plus_h": self.m_plus_h,
            "m_plus_na": self.m_plus_na,
            # "origin_reference": self.origin_reference,
            # "npclassifier": self.npclassifier,
            # "external_ids": self.external_ids,
            "ncbi_taxid": self.ncbi_taxid,
            "genus": self.genus,
            "species": self.species,
            # "gnps": self.gnps,
            # "npmrd": self.npmrd,
            # "mibig": self.mibig,
            "classyfire_class": self.classyfire_class,
            "classyfire_subclass": self.classyfire_subclass,
        }





class NPAtlas(ExternalBaseClass):

    def __init__(self, atlas_json_path) -> None:
        super().__init__()
        self.path = atlas_json_path
        self.entries = []


    def _hydrate(self):
        with open(self.path, "r") as f:
            for i in json.load(f):
                self.entries.append(NPAtlasEntry(i))

    def add_nodes_to_neo4j(self):
        self._add_to_neo4j(
            """
                CREATE CONSTRAINT npatlas_entry IF NOT EXISTS
                FOR (n:npatlas)
                REQUIRE (n.uid) IS UNIQUE;
            """
        )
        for i in batched((i._node_prop_dict for i in self.entries), 10000):
            self._add_to_neo4j(
                """
                WITH $input as input
                UNWIND input as row
                MERGE (np:npatlas {uid: row.uid})
                ON CREATE SET np+= {
                    original_name: row.original_name,
                    mol_formula: row.mol_formula,
                    mol_weight: row.mol_weight,
                    exact_mass: row.exact_mass,
                    inchikey: row.inchikey,
                    smiles: row.smiles,
                    cluster_id: row.cluster_id,
                    node_id: row.node_id,
                    synonyms: row.synonyms,
                    inchi: row.inchi,
                    m_plus_h: row.m_plus_h,
                    m_plus_na: row.m_plus_na,
                    //origin_reference: row.origin_reference,
                    //npclassifier: row.npclassifier,
                    //external_ids: row.external_ids,
                    ncbi_taxid: row.ncbi_taxid,
                    genus: row.genus,
                    species: row.species,
                    //gnps: row.gnps,
                    //npmrd: row.npmrd,
                    //mibig: row.mibig,
                    classyfire_class: row.classyfire_class,
                    classyfire_subclass: row.classyfire_subclass
                }
                """,
                input=i,
            )




    def merge_with_mibig(self):
        for mibig_obj in self.mibig:
            self._add_to_neo4j(
                """
                MERGE (np:npatlas {uid: $uid})
                MERGE (bgc:assembly {uid: $mibig})
                MERGE (bgc)-[:PRODUCES]->(np)
                """,
                uid=self.uid,
                mibig=mibig_obj.uid,
            )

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

    def connect_to_chem(self):
        ChemicalCompound()
