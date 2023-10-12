from typing import Dict, Set

from socialgene.hashing.hashing import hasher
from socialgene.neo4j.neo4j import GraphDriver
from socialgene.utils.logging import log


class Moiety:
    __slots__ = ["name"]

    def __init__(self, name) -> None:
        self.name = name.lower()


class Target:
    __slots__ = ["name"]

    def __init__(self, name) -> None:
        self.name = name


class Activity:
    __slots__ = ["name"]

    def __init__(self, name) -> None:
        self.name = name.lower()


class Compound:
    __slots__ = [
        "smiles",
        "name",
        "xrefs",
        "mol_mass",
        "mol_formula",
        "moieties",
        "targets",
        "activities",
        "uid",
        "source",
    ]

    def __init__(
        self,
        smiles: str = None,
        name: str = None,
        xrefs: Dict = {},
        mol_mass: float = None,
        mol_formula: str = None,
        moieties: Set = None,
        targets: Set = None,
        activities: Set = None,
        uid: str = None,
        source: str = "unknown",
    ):
        self.smiles = smiles
        self.name = name
        self.xrefs = xrefs
        self.mol_mass = mol_mass
        self.mol_formula = mol_formula
        self.moieties = moieties
        self.targets = targets
        self.activities = activities
        self.uid = uid
        self.source = source

    def write_all(self):
        self._db_merge_chem()
        self._db_link_xrefs()
        self._db_link_targets()
        self._db_link_moieties()
        self._db_link_activities()

    def _db_merge_chem(self, source: str = None):
        if self.uid:
            # TODO: this will overwrite the properties of an existing node
            # But can't merge with everything because a single prop difference
            # will create an entirely new node
            self._add_to_neo4j(
                """
                MERGE (chem:chemical {uid: $uid})
                SET chem = $props
                """,
                uid=self.uid,
                props={
                    "uid": self.uid,
                    "mol_mass": self.mol_mass,
                    "mol_formula": self.mol_formula,
                    "name": self.name,
                    "smiles": self.smiles,
                },
            )

    def _db_link_xrefs(self):
        # TODO: this should rely on the other db methods (currently not implemented)
        # e.g. npatlas ids should use npatlas class
        if self.uid and self.xrefs:
            for xref_key, xref_value in self.xrefs.items():
                # can't  pass 'xref_key' as param due to limitation of cypher
                self._add_to_neo4j(
                    f"""
                    MERGE (xref:{xref_key} {{uid: $input_value}})
                    MERGE (chem:chemical {{uid: $uid}})
                    MERGE (chem)-[:CROSSREF {{source:$source}}]->(xref)
                    """,
                    input_value=str(xref_value),
                    uid=self.uid,
                    source=self.source,
                )

    def _db_link_targets(self):
        if self.uid and self.targets:
            self._add_to_neo4j(
                """
                    WITH $targets as targets
                    UNWIND targets AS target
                    MERGE (t1:target {uid: target})
                    MERGE (chem:chemical {uid: $uid})
                    MERGE (chem)-[:ACTS_ON {source:$source}]->(t1)
                    """,
                targets=[i.name for i in self.targets],
                uid=self.uid,
                source=self.source,
            )

    def _db_link_moieties(self):
        if self.uid and self.moieties:
            self._add_to_neo4j(
                """
                    WITH $moieties as moieties
                    UNWIND moieties AS moiety
                    MERGE (t1:moiety {uid: moiety})
                    MERGE (chem:chemical {uid: $uid})
                    MERGE (chem)<-[:SUBSTRUCTURE {source:$source}]->(t1)
                    """,
                moieties=[i.name for i in self.moieties],
                uid=self.uid,
                source=self.source,
            )

    def _db_link_activities(self):
        if self.uid and self.activities:
            self._add_to_neo4j(
                """
                    WITH $activities as activities
                    UNWIND activities AS activity
                    MERGE (t1:activity {uid: activity})
                    MERGE (chem:chemical {uid: $uid})
                    MERGE (chem)-[:SHOWS {source:$source}]->(t1)
                    """,
                activities=[i.name for i in self.activities],
                uid=self.uid,
                source=self.source,
            )

    def _add_to_neo4j(self, statement, **kwargs):
        summary = (
            GraphDriver()
            .driver.execute_query(
                statement,
                **kwargs,
                database_="neo4j",
            )
            .summary
        )

        if summary.metadata.get("stats"):
            log.info(
                f"{summary.metadata.get('stats').get('properties-set')} properties modified"
            )
        else:
            log.info("No properties modified")


class Mibig_Compound(Compound):
    def __init__(self, cmpd_dict: Dict):
        """Create a Compound class using extracted "compounds" json dict from mibig json

        Args:
            cmpd_dict (Dict): extracted dict  `for cmpd_dict in json.load(mibig_bgc_json_file)["cluster"]["compounds"]:`
        """
        super().__init__()
        self.cmpd_dict = cmpd_dict

    def assign(self):
        self._assign_smiles()
        self._assign_name()
        self._assign_xref()
        self._assign_mol_mass()
        self._assign_mol_formula()
        self._assign_moieties()
        self._assign_targets()
        self._assign_activities()
        self._assign_uid()
        self.source = "mibig"

    def _assign_uid(self):
        # This still works even if None; note: CRC64 of double None is 6B2273F47D2273F4
        self.uid = hasher(f"{self.name.lower()}{self.smiles}")

    def _assign_activities(self):
        if "chem_acts" in self.cmpd_dict:
            self.activities = {
                Activity(name=i["activity"]) for i in self.cmpd_dict.get("chem_acts")
            }

    def _assign_moieties(self):
        if "chem_moieties" in self.cmpd_dict:
            self.moieties = {
                Moiety(name=i["moiety"].lower())
                for i in self.cmpd_dict.get("chem_moieties")
            }

    def _assign_smiles(self):
        if "chem_struct" in self.cmpd_dict:
            self.smiles = self.cmpd_dict["chem_struct"]

    def _assign_targets(self):
        if "chem_targets" in self.cmpd_dict:
            self.targets = {
                Target(name=i["target"])
                for i in self.cmpd_dict.get("chem_targets")
                if "target" in i
            }

    def _assign_name(self):
        if "compound" in self.cmpd_dict:
            self.name = self.cmpd_dict.get("compound")

    def _assign_xref(self):
        # gets db by name and puts into dict
        # e.g. "npatlas:NPA01961", becomes {"npatlas": "NPA01961"}
        if "database_id" in self.cmpd_dict:
            for xref in self.cmpd_dict.get("database_id"):
                temp = xref.split(":")
                self.xrefs = self.xrefs | {temp[0]: temp[1]}

    def _assign_mol_mass(self):
        if "mol_mass" in self.cmpd_dict:
            self.mol_mass = self.cmpd_dict.get("mol_mass")

    def _assign_mol_formula(self):
        if "molecular_formula" in self.cmpd_dict:
            self.mol_formula = self.cmpd_dict.get("molecular_formula")
