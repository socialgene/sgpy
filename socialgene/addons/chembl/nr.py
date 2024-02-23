"""https://www.ebi.ac.uk/chembl"""

from socialgene.neo4j.neo4j_element import Node


class CHEMBL(Node):
    neo4j_label = ["chembl"]
    description = "Represents a CHEMBL term"
    property_specification = {
        "uid": str,
    }
    constraints_unique = ["uid"]
