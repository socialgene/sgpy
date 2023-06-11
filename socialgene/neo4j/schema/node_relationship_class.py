from typing import List

from abc import ABC, abstractmethod
from collections.abc import Generator


class Neo4jElement:
    def __init__(
        self,
        neo4j_label: str,
        header_filename: str,
        target_subdirectory: str,
        target_extension: str,
        header: List[str],
        multilabel: bool = False,
    ):
        """This defines nodes and relationships that will be imported  into Neo4j  with the Neo4j Admin Import tool
        For more information, especially about what makes up the headers, see: https://neo4j.com/docs/operations-manual/current/tutorial/neo4j-admin-import/

        Args:
            neo4j_label (str): this will become the Neo4j node or relationship LABEL
            header_filename (str): the name of the header file used for Neo4j admin import (basically column names, but not quite)
            target_subdirectory (str): subdirectory the import data can be found in  (e.g. for non-redundant protein nodes it would be 'protein_info' because data is within: `$outdir/socialgene_neo4j/import/protein_info`)
            target_extension (str): extension that is unique to the wanted data files (scoped within target_subdirectory)
            header (List): list of strings, each string will make up a column in a tab separated header file for Neo4j admin import (basically column names, but not quite)
            multilabel(bool): are LABELS specified within the tsv?
        """
        self.neo4j_label = neo4j_label
        self.header_filename = header_filename
        self.target_subdirectory = target_subdirectory
        self.target_extension = target_extension
        self.header = header
        self.multilabel = multilabel

    def __hash__(self):
        """Node or relationship data is defined as the unique combination of:
        neo4j LABEL and the file that LABEL data will be imported from
        """
        return hash((self.neo4j_label, self.target_subdirectory, self.target_extension))

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return NotImplemented
        return (
            self.neo4j_label == other.neo4j_label
            and self.target_subdirectory == other.target_subdirectory
            and self.target_extension == other.target_extension
        )


class NR(ABC):
    def __init__(
        self,
    ):
        # For "set()" info see the "__hash__" function in the Neo4jElement() class
        self.nodes = set()
        self.relationships = set()
        self._import()

    def add_node(self, **kwargs):
        self.nodes.add(Neo4jElement(**kwargs))

    def add_relationship(self, **kwargs):
        self.relationships.add(Neo4jElement(**kwargs))

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
