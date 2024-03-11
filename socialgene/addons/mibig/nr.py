"""https://mibig.secondarymetabolites.org/"""

import re

from socialgene.addons.chemistry.nr import ChemicalCompoundNode
from socialgene.neo4j.neo4j_element import Node, Relationship
from socialgene.nextflow.nodes import ASSEMBLY


class Mibig(ASSEMBLY):
    neo4j_label = ["assembly", "mibig"]
    description = "Represents a single Mibig entry"
    property_specification = {
        "uid": str,
    }
    required_properties = ["uid"]

    def __init__(self, properties=None) -> None:
        super().__init__()
        if properties:
            self.properties["uid"] = self.attempt_to_fix_uid(properties["uid"])

    @staticmethod
    def attempt_to_fix_uid(uid):
        # npatlas has some mibig ids prefixed, some not, try to add prefix here
        # https://github.com/NPAtlas/npatlas_website/issues/13
        if uid.startswith("BGC"):
            temp = uid
        else:
            temp = "BGC" + uid.zfill(7)
        # Validate uid
        if re.match("BGC[0-9]{7}", temp):
            return temp
        else:
            raise ValueError(f"Unexpected mibig ID {uid}")


class Substrate(ChemicalCompoundNode):
    neo4j_label = ChemicalCompoundNode.neo4j_label + ["substrate"]
    description = "Mibig substrate (e.g. NRPS monomer)"


class MibigCompound(Node):
    neo4j_label = ["mibig_compound"]
    description = "Represents a single Mibig compound"
    property_specification = {
        "name": str,
        "smiles": str,
        "inchi": str,
    }
    required_properties = ["uid"]
    constraints_unique = ["uid"]


class MibigActivity(Node):
    neo4j_label = ["mibig_activity"]
    description = "Represents a single Mibig bioactivity"
    property_specification = {
        "uid": str,
    }
    required_properties = ["uid"]
    constraints_unique = ["uid"]


class Mibig_Biosynthetic_Class(Node):
    neo4j_label = ["mibig_biosynthetic_class"]
    description = "Represents a single Mibig biosynthetic class"
    property_specification = {
        "uid": str,
    }
    required_properties = ["uid"]
    constraints_unique = ["uid"]


class MibigToMibigCompound(Relationship):
    neo4j_label = "PRODUCES"
    description = "Connects a Mibig BGC to a Mibig compound"
    start_class = Mibig
    end_class = MibigCompound
