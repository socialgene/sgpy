"""https://www.ebi.ac.uk/chebi/"""

from socialgene.neo4j.neo4j_element import Node


class ChebiNode(Node):
    def __init__(self, *args, **kwargs):
        super().__init__(
            neo4j_label="chebi",
            description="Represents a ChEBI term",
            properties={
                "uid": str,
            },
        )
