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
            self.__properties = header


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

    def __rich_console__(
        self, console: Console, options: ConsoleOptions
    ) -> RenderResult:  # pragma: no cover
        table = Table(title="Node")
        table.add_column("Label", justify="left", style="cyan", no_wrap=True, ratio=1)
        table.add_column(
            "Description",
            justify="left",
            style="cyan",
            no_wrap=False,
            ratio=4,
            max_width=50,
        )
        table.add_column("Nextflow results subdirectory", style="magenta", ratio=1)
        table.add_column("Neo4j header file", style="magenta", ratio=1)
        table.add_column("sdsdsd", style="magenta", ratio=1)
        table.add_row(
            self.neo4j_label,
            self.target_subdirectory,
            self.target_subdirectory,
            self.header_filename,
            self.header,
        )
        yield table


class Relationship(Neo4jElement):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

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
        return f"(:{self._start})-[:{self._Neo4jElement__neo4j_label}]->(:{self._end})"
