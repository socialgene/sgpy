"""http://classyfire.wishartlab.com"""

from socialgene.addons.chebi.nr import ChebiNode
from socialgene.neo4j.neo4j_element import Node, Relationship


class ClassyFireNode(Node):
    neo4j_label = ["classyfire"]
    description = "Represents a classyfire chemical ontology term"
    property_specification = {"uid": int, "name": str, "definition": str}
    required_properties = ["uid"]
    constraints_unique = ["uid"]


class ClassyFireIsA(Relationship):
    neo4j_label = "IS_A"
    description = "Represents a relationship between chemont nodes"
    start_class = ClassyFireNode
    end_class = ClassyFireNode


class ClassyFireSynonym(Relationship):
    neo4j_label = "SYNONYM"
    description = "Represents a synonym relationship between chemont chebi"
    start_class = ClassyFireNode
    end_class = ChebiNode
