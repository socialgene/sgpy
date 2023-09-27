from abc import ABC, abstractmethod
from collections.abc import Generator

from socialgene.neo4j.schema.neo4j_element import Neo4jElement


class NR(ABC):
    def __init__(
        self,
    ):
        # For "set()" info see the "__hash__" function in the Neo4jElement() class
        self.nodes = set()
        self.relationships = set()
        self._import()

    # def add_node(self, **kwargs):
    #     self.nodes.add(Neo4jElement(**kwargs))

    # # def add_relationship(self, **kwargs):
    #     self.relationships.add(Neo4jElement(**kwargs))

    def _get_by_label(self, x, y):
        return (i for i in x if i.neo4j_label in y)

    def get_nodes(self, input):
        return self._get_by_label(x=self.nodes, y=input)

    def get_relationships(self, input):
        return self._get_by_label(x=self.relationships, y=input)

    def get_relationship_by_neo4j_label(
        self, neo4j_label: str
    ) -> Generator[Neo4jElement]:
        """Filter relationships by neo4j label and return generator

        Args:
            neo4j_label (str): neo4j_label

        Yields:
            Generator[Neo4jElement]: Neo4jElement relationships that matched input label
        """
        return (i for i in self.relationship if i.neo4j_label == neo4j_label)

    @abstractmethod
    def _import():
        pass
