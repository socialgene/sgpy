from typing import List
from socialgene.utils.logging import log


class Neo4jElement:
    def __init__(
        self,
        neo4j_label: str,
        header_filename: str,
        target_subdirectory: str,
        target_extension: str,
        header: List[str],
    ):
        """This defines nodes and relationships that will be imported  into Neo4j  with the Neo4j Admin Import tool
        For more information, especially about what makes up the headers, see: https://neo4j.com/docs/operations-manual/current/tutorial/neo4j-admin-import/

        Args:
            neo4j_label (str): this will become the Neo4j node or relationshipLABEL
            header_filename (str): the name of the header file used for Neo4j admin import (basically column names, but not quite)
            target_subdirectory (str): subdirectory the import data can be found in  (e.g. for non-redundant protein nodes it would be 'protein_info' because data is within: `$outdir/socialgene_neo4j/import/protein_info`)
            target_extension (str): extension that is unique to the wanted data files (scoped within target_subdirectory)
            header (List): list of strings, each string will make up a column in a tab separated header file for Neo4j admin import (basically column names, but not quite)
        """
        self.neo4j_label = neo4j_label
        self.header_filename = header_filename
        self.target_subdirectory = target_subdirectory
        self.target_extension = target_extension
        self.header = header

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
