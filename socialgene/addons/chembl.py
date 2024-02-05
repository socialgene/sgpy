"""https://www.ebi.ac.uk/chembl"""

from socialgene.neo4j.neo4j_element import Node


class CHEMBL(Node):
    def __init__(self, *args, **kwargs):
        super().__init__(
            neo4j_label="chembl",
            description="Represents a CHEMBL term",
            properties={
                "uid": "string",
            },
        )
