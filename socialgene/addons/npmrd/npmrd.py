"""https://np-mrd.org/"""


class Npmrd:
    neo4j_label = "npmrd"
    description = "Represents a single NP-MRD entry"
    property_specification = {
        "uid": str,
    }
    required_properties = ["uid"]
