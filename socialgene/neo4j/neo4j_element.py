from abc import ABC
from typing import List
import re
from rich.console import Console, ConsoleOptions, RenderResult
from rich.table import Table
from socialgene.utils.logging import log


class Neo4jElement(ABC):
    def __init__(
        self,
        neo4j_label: str = None,
        description: str = None,
        header_filename: str = None,
        target_subdirectory: str = None,
        target_extension: str = None,
        multilabel: bool = False,
        header: List[str] = None,
        properties: List[str] = None,
        **kwargs,
    ):
        """This defines nodes and relationships that will be imported  into Neo4j  with the Neo4j Admin Import tool
        For more information, especially about what makes up the headers, see: https://neo4j.com/docs/operations-manual/current/tutorial/neo4j-admin-import/

        Args:
            neo4j_label (str): this will become the Neo4j node or relationship LABEL
            description (str): helpful information about what the element represents
            header_filename (str): the name of the header file used for Neo4j admin import (basically column names, but not quite)
            target_subdirectory (str): subdirectory the import data can be found in  (e.g. for non-redundant protein nodes it would be 'protein_info' because data is within: `$outdir/socialgene_neo4j/import/protein_info`)
            target_extension (str): extension that is unique to the wanted data files (scoped within target_subdirectory)
            header (List): list of strings, each string will make up a column in a tab separated header file for Neo4j admin import (basically column names, but not quite)
            multilabel(bool): are LABELS specified within the tsv?
        """
        # TODO: multi-label
        self.__neo4j_label = neo4j_label
        self.__description = description
        self.__header_filename = header_filename
        self.__target_subdirectory = target_subdirectory
        self.__target_extension = target_extension
        self.__header = header
        self.__properties = properties
        self.__multilabel = multilabel
        if self.__properties is None:
            if self.__header is not None:
                self.__properties = {i: "string" for i in header}

    def __hash__(self):
        return hash((self.__neo4j_label))


class Node(Neo4jElement):
    """Represents a single Node"""

    def __init__(
        self,
        contraints: List[str] = None,
        contraints_unique: List[str] = None,
        *args,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.__contraints = contraints
        self.__contraints_unique = contraints_unique

    def __hash__(self):
        return hash((self.neo4j_label))


class Relationship(Neo4jElement):
    def __init__(self, start=None, end=None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.start = start
        self.end = end
        if not self.start and not self.end:
            if self._Neo4jElement__header is not None:
                self.start = self._start
                self.end = self._end

    def _start_field(self):
        for i in self._Neo4jElement__header:
            if "START_ID" in i:
                return i

    def _end_field(self):
        for i in self._Neo4jElement__header:
            if "END_ID" in i:
                return i

    @property
    def _start(self):
        try:
            return re.findall(r"\((.*?)\)", self._start_field())[0]
        except Exception as e:
            log.debug(e)

    @property
    def _end(self):
        try:
            return re.findall(r"\((.*?)\)", self._end_field())[0]
        except Exception as e:
            log.debug(e)

    @property
    def _cypher_string(self):
        if self.start and self.end:
            return f"(:{self.start()._Neo4jElement__neo4j_label})-[:{self._Neo4jElement__neo4j_label}]->(:{self.end()._Neo4jElement__neo4j_label})"
