"""http://classyfire.wishartlab.com"""

from socialgene.addons.chebi.nr import ChebiNode
from socialgene.neo4j.neo4j_element import Node, Relationship
from socialgene.utils.logging import log


class ClassyFireNode(Node):
    neo4j_label = "classyfire"
    description = "Represents a CHEMBL term"
    property_specification = {
        "uid": str,
    }

class ClassyChemontNode(Node):
    neo4j_label = "chemont"
    description = "Represents a classyfire chemical ontology term"
    property_specification = {"uid": str, "name": str, "definition": str}


class ClassyFireSynonym(Relationship):
    neo4j_label = "SYNONYM"
    description = "Represents a synonym relationship between chemont chebi"
    start_class = ClassyChemontNode
    end_class = ChebiNode


class ClassyFireIsA(Relationship):
    neo4j_label = "IS_A"
    description = "Represents a relationship between chemont nodes"
    start_class = ClassyChemontNode
    end_class = ClassyChemontNode

