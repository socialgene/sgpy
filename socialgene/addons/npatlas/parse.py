from pathlib import Path

from socialgene.addons.chebi.nr import ChebiNode
from socialgene.addons.chemistry.nr import ChemicalCompoundNode
from socialgene.addons.classyfire.nr import ClassyFireNode
from socialgene.addons.gnps_library.nr import GnpsLibrarySpectrumNode
from socialgene.addons.mibig.nr import Mibig
from socialgene.addons.npatlas.nr import (
    MibigToNPAtlas,
    NPAtlasNode,
    NPAtlasPublication,
    NPAtlasToChebi,
    NPAtlasToChem,
    NPAtlasToClassyFireAlternativeParents,
    NPAtlasToClassyFireDirectParent,
    NPAtlasToClassyFireIntermediateNodes,
    NPAtlasToClassyFireLowestClass,
    NPAtlasToGnps,
    NPAtlasToNpclassifierClass,
    NPAtlasToNpclassifierPathway,
    NPAtlasToNpclassifierSuperclass,
    NPAtlasToPublication,
    NPAtlasToTaxID,
)
from socialgene.addons.npclassifier.nr import (
    NPClassifierClass,
    NPClassifierPathway,
    NPClassifierSuperclass,
)
from socialgene.addons.npmrd.nr import Npmrd
from socialgene.base.chem import ChemicalCompound
from socialgene.nextflow.nodes import TAXID
from socialgene.utils.download import download as downloader
from socialgene.utils.logging import log

NPATALAS_URL = "https://www.npatlas.org/static/downloads/NPAtlas_download.json"


def _download_npatlas(outpath, url=NPATALAS_URL):
    # check if exist and not empty
    if Path(outpath).exists() and Path(outpath).stat().st_size > 0:
        log.debug(f"File already exists at {outpath}")
    else:
        downloader(url, outpath)


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
        "external_ids",
        "ncbi_taxid",
        "genus",
        "species",
        "lowest_classyfire",
        "gnps_ids",
        "mibig_ids",
        "npmrd_ids",
        "predicted_chebi_terms",
        "classyfire_direct_parent",
        "classyfire_intermediate_nodes",
        "classyfire_alternative_parents",
        "npclassifier_classes",
        "npclassifier_pathways",
        "npclassifier_superclasses",
        "npclassifier_glycoside",
    ]

    def __init__(self, entry) -> None:
        self.entry = entry
        self.uid = None
        self.original_name = None
        self.mol_formula = None
        self.mol_weight = None
        self.exact_mass = None
        self.inchikey = None
        self.smiles = None
        self.cluster_id = None
        self.node_id = None
        self.synonyms = None
        self.inchi = None
        self.m_plus_h = None
        self.m_plus_na = None
        self.origin_reference = None
        self.external_ids = None
        self.ncbi_taxid = None
        self.genus = None
        self.species = None
        self.lowest_classyfire = None
        self.gnps_ids = set()
        self.mibig_ids = set()
        self.npmrd_ids = set()
        self.predicted_chebi_terms = {}
        self.classyfire_direct_parent = None
        self.classyfire_intermediate_nodes = None
        self.classyfire_alternative_parents = None
        self.npclassifier_classes = None
        self.npclassifier_pathways = None
        self.npclassifier_superclasses = None
        self.npclassifier_glycoside = None

    def parse(self):
        self._parse_single_entry()
        self._assign_lowest_classyfire()
        self._assign_classyfire_direct_parent()
        self._assign_classyfire_intermediate_nodes()
        self._assign_classyfire_alternative_parents()
        self._assign_npclassifier()
        self._assign_publication()

    def _parse_single_entry(self):
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
        self.origin_reference = self.entry.get("origin_reference", None)
        self.synonyms = None
        self._assign_taxon()
        self._assign_external_ids()
        self._assign_chebi_terms()

    def _assign_publication(self):
        self.origin_reference = NPAtlasPublication()
        try:
            self.origin_reference.properties["doi"] = (
                self.origin_reference._extract_doi(
                    self.entry.get("origin_reference").get("doi", None)
                )
            )
        except Exception as e:
            log.debug(f"Failed to extract doi: {e}")

        for i in ["pmid", "authors", "title", "journal", "year"]:
            try:
                self.origin_reference.properties[i] = self.entry.get(
                    "origin_reference"
                ).get(i, None)
            except Exception as e:
                log.debug(f"Failed to extract {i}: {e}")

    def _assign_lowest_classyfire(self):
        # not all entries have classyfire data at every level, just assign the most specific one
        for i in ["kingdom", "superclass", "class", "subclass"]:
            try:
                self.lowest_classyfire = self.entry["classyfire"][i][
                    "chemont_id"
                ].removeprefix("CHEMONTID:")
            except Exception as e:
                log.debug(e)

    def _assign_classyfire_direct_parent(self):
        try:
            self.classyfire_direct_parent = self.entry["classyfire"]["direct_parent"][
                "chemont_id"
            ].removeprefix("CHEMONTID:")
            self.classyfire_direct_parent = int(self.classyfire_direct_parent)
        except Exception as e:
            log.debug(e)

    def _assign_classyfire_intermediate_nodes(self):
        try:
            self.classyfire_intermediate_nodes = [
                int(x["chemont_id"].removeprefix("CHEMONTID:"))
                for x in self.entry["classyfire"]["intermediate_nodes"]
            ]
        except Exception as e:
            log.debug(e)

    def _assign_classyfire_alternative_parents(self):
        try:
            self.classyfire_alternative_parents = [
                int(x["chemont_id"].removeprefix("CHEMONTID:"))
                for x in self.entry["classyfire"]["alternative_parents"]
            ]
        except Exception as e:
            log.debug(e)

    def _assign_npclassifier(self):
        try:
            self.npclassifier_classes = self.entry["npclassifier"]["class_results"]
        except Exception as e:
            log.debug(e)
        try:
            self.npclassifier_pathways = self.entry["npclassifier"]["pathway_results"]
        except Exception as e:
            log.debug(e)
        try:
            self.npclassifier_superclasses = self.entry["npclassifier"][
                "superclass_results"
            ]
        except Exception as e:
            log.debug(e)
        try:
            if self.entry["npclassifier"]["isglycoside"]:
                self.npclassifier_glycoside = True
        except Exception as e:
            log.debug(e)

    def _assign_taxon(self):
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

    def _assign_chebi_terms(self):
        try:
            for term in self.entry["classyfire"]["predicted_chebi_terms"]:
                # extract name and uid as int from strings like "hydroxycoumarin (CHEBI:37912)"
                uid = term.split(" ")[-1].removeprefix("(CHEBI:").removesuffix(")")
                name = term.split("(")[0]
                uid = uid.strip()
                name = name.strip()
                self.predicted_chebi_terms[uid] = name
        except Exception as e:
            log.debug(e)

    def _assign_external_ids(self):
        for external_id in self.entry["external_ids"]:
            try:
                match external_id["external_db_name"]:
                    case "gnps":
                        # get first becasue some npaid have multiple gnps ids, eg:
                        # "CCMSLIB00005722620%NCGC00179860-02!6-[(3R,4S,5S,7R)-7-[(2S,3S,5S)-5-ethyl-5-[(2R,5R,6S)-5-ethyl-5-hydroxy-6-methyloxan-2-yl]-3-methyloxolan-2-yl]-4-hydroxy-3,5-dimethyl-6-oxononyl]-2-hydroxy-3-methylbenzoic acid [IIN-based on: CCMSLIB00000848016]%3""
                        uid = GnpsLibrarySpectrumNode._extract_all_CCMSLIB(
                            external_id["external_db_code"]
                        )[0]
                        self.gnps_ids.add(
                            GnpsLibrarySpectrumNode(properties={"uid": uid})
                        )
                    case "npmrd":
                        self.npmrd.add(
                            Npmrd(properties={"uid": external_id["external_db_code"]})
                        )
                    case "mibig":
                        self.mibig_ids.add(
                            Mibig(properties={"uid": external_id["external_db_code"]})
                        )
                    case _:
                        pass
            except Exception as e:
                log.debug(e)

    @property
    def node(self):
        return NPAtlasNode(
            properties={
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
                "genus": self.genus,
                "species": self.species,
            }
        )

    def link_publication(self) -> set:
        if self.origin_reference is None:
            return set()
        return {NPAtlasToPublication(start=self.node, end=self.origin_reference)}

    def link_ncbi_taxid(self) -> set:
        if self.ncbi_taxid is None:
            return set()
        return {
            NPAtlasToTaxID(
                start=TAXID(properties={"uid": str(self.ncbi_taxid)}), end=self.node
            )
        }

    def link_gnps(self) -> set:
        if len(self.gnps_ids) == 0:
            return set()
        return {NPAtlasToGnps(start=self.node, end=x) for x in self.gnps_ids}

    def link_mibig(self) -> set:
        if len(self.mibig_ids) == 0:
            return set()
        return {MibigToNPAtlas(start=x, end=self.node) for x in self.mibig_ids}

    def link_classyfire(self) -> set:
        if self.lowest_classyfire is None:
            return set()
        return {
            NPAtlasToClassyFireLowestClass(
                start=self.node,
                end=ClassyFireNode(properties={"uid": int(self.lowest_classyfire)}),
            )
        }

    def link_classyfire_direct_parent(self) -> set:
        if self.classyfire_direct_parent is None:
            return set()
        return {
            NPAtlasToClassyFireDirectParent(
                start=self.node,
                end=ClassyFireNode(
                    properties={"uid": int(self.classyfire_direct_parent)}
                ),
            )
        }

    def link_classyfire_intermediate_nodes(self) -> set:
        if self.classyfire_intermediate_nodes is None:
            return set()
        return {
            NPAtlasToClassyFireIntermediateNodes(
                start=self.node, end=ClassyFireNode(properties={"uid": int(x)})
            )
            for x in self.classyfire_intermediate_nodes
        }

    def link_classyfire_alternative_parents(self) -> set:
        if self.classyfire_alternative_parents is None:
            return set()
        return {
            NPAtlasToClassyFireAlternativeParents(
                start=self.node, end=ClassyFireNode(properties={"uid": int(x)})
            )
            for x in self.classyfire_alternative_parents
        }

    def link_npclassifier_classes(self) -> set:
        if self.npclassifier_classes is None:
            return set()
        return {
            NPAtlasToNpclassifierClass(
                start=self.node, end=NPClassifierClass(properties={"uid": x})
            )
            for x in self.npclassifier_classes
        }

    def link_npclassifier_pathways(self) -> set:
        if self.npclassifier_pathways is None:
            return set()
        return {
            NPAtlasToNpclassifierPathway(
                start=self.node, end=NPClassifierPathway(properties={"uid": x})
            )
            for x in self.npclassifier_pathways
        }

    def link_npclassifier_superclasses(self) -> set:
        if self.npclassifier_superclasses is None:
            return set()
        return {
            NPAtlasToNpclassifierSuperclass(
                start=self.node, end=NPClassifierSuperclass(properties={"uid": x})
            )
            for x in self.npclassifier_superclasses
        }

    def link_chem(self) -> set:
        try:
            chem = None
            if self.inchi:
                chem = ChemicalCompound(self.inchi)
            elif self.smiles:
                chem = ChemicalCompound(self.smiles)
            if chem:
                a = ChemicalCompoundNode()
                a.fill_from_dict(chem.base_properties | chem.hash_dict)
                return {NPAtlasToChem(self.node, a)}
        except Exception:
            return set()

    def link_chebi(self) -> set:
        try:
            return {
                NPAtlasToChebi(
                    start=self.node,
                    end=ChebiNode(properties={"uid": int(k), "name": v}),
                )
                for k, v in self.predicted_chebi_terms.items()
            }
        except Exception:
            return set()

    def get_links(self):
        return {
            NPAtlasToPublication: self.link_publication(),
            NPAtlasToTaxID: self.link_ncbi_taxid(),
            NPAtlasToGnps: self.link_gnps(),
            MibigToNPAtlas: self.link_mibig(),
            NPAtlasToClassyFireLowestClass: self.link_classyfire(),
            NPAtlasToClassyFireDirectParent: self.link_classyfire_direct_parent(),
            NPAtlasToClassyFireIntermediateNodes: self.link_classyfire_intermediate_nodes(),
            NPAtlasToClassyFireAlternativeParents: self.link_classyfire_alternative_parents(),
            NPAtlasToNpclassifierClass: self.link_npclassifier_classes(),
            NPAtlasToNpclassifierPathway: self.link_npclassifier_pathways(),
            NPAtlasToNpclassifierSuperclass: self.link_npclassifier_superclasses(),
            NPAtlasToChem: self.link_chem(),
            NPAtlasToChebi: self.link_chebi(),
        }
