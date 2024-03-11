"""https://np-mrd.org/"""

from socialgene.neo4j.neo4j_element import Node


class Npmrd(Node):
    neo4j_label = ["npmrd"]
    description = "Represents a single NP-MRD entry"
    property_specification = {
        "uid": str,
    }
    required_properties = ["uid"]
    constraints_unique = ["uid"]
