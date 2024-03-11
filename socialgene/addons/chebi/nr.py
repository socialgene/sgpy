"""https://www.ebi.ac.uk/chebi/"""

from socialgene.neo4j.neo4j_element import Node


class ChebiNode(Node):
    neo4j_label = ["chebi"]
    description = "Represents a ChEBI term"
    property_specification = {
        "uid": int,
        "name": str,
    }
    constraints_unique = ["uid"]
