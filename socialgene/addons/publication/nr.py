import re

from socialgene.neo4j.neo4j_element import Node


class Publication(Node):
    neo4j_label = ["publication"]
    description = "Represents a publication"
    property_specification = {
        "doi": str,
        "pmid": int,
        "authors": str,
        "title": str,
        "journal": str,
        "year": int,
    }
    required_properties = ["doi"]
    constraints_unique = ["doi"]

    @staticmethod
    def _extract_doi(input):
        return re.search("10.\\d{4,9}/[-._;()/:a-z0-9A-Z]+", input).group()
