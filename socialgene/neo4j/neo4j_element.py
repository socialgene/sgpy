from abc import ABC, abstractmethod
from typing import List
import re
from rich.console import Console, ConsoleOptions, RenderResult
from rich.table import Table
from socialgene.neo4j.neo4j import GraphDriver
from socialgene.utils.logging import log
from itertools import batched


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
        property_specification: dict = None,
        required_properties: List[str] = None,
        properties: dict = None,
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
        self.__neo4j_label = neo4j_label
        self.__description = description
        self.__header_filename = header_filename
        self.__target_subdirectory = target_subdirectory
        self.__target_extension = target_extension
        self.__header = header
        self.__property_specification = property_specification
        if self.__property_specification is None:
            if self.__header is not None:
                self.__property_specification = {i: "string" for i in header}
        # properties aren't actually required
        if self.__property_specification is None:
            self.__property_specification = {}
        # if not provided then all properties are required
        self.__required_properties = required_properties
        if self.__required_properties is None:
            self.__required_properties = list(self.__property_specification.keys())
        if self.__required_properties is None:
            self.__required_properties = []
        self.__multilabel = multilabel
        if properties is None:
            self.__properties = {}
        if properties is not None:
            self.__properties = properties
            self._clean_properties()


    @property
    def __required_properties_dict(self):
        return {i: self.__properties.get(i) for i in self.__required_properties}

    @property
    def __optional_properties_dict(self):
        return {
            i: self.__properties.get(i)
            for i in self.__properties
            if i not in self.__required_properties
        }

    def __hash__(self):
        return hash((self.__neo4j_label))

    def _check_required_properties(
        self,
    ):
        missing = [i for i in self.__required_properties if i not in self.__properties]
        if missing:
            raise ValueError(f"{self.__class__} missing properties: {missing}")
        return missing

    def _check_property_types(self):
        bad = {
            k: f"expected {self.__property_specification.get(k).__name__} but got {type(v).__name__}"
            for k, v in self.__properties.items()
            if not isinstance(v, self.__property_specification.get(k))
        }
        if bad:
            raise ValueError(f"{self.__class__} bad properties:\n\t\t{bad}")
        self.__properties = {
            k: v
            for k, v in self.__properties.items()
            if isinstance(v, self.__property_specification.get(k))
        }

    def _clean_properties(
        self,
    ):
        self._check_required_properties()
        self._check_property_types()


class Node(Neo4jElement):
    """Represents a single Node"""

    def __init__(
        self,
        uid: List[str] = None,
        constraints: List[str] = None,
        unique_constraints: List[str] = None,
        nonunique_index: List[str] = None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.__uid = uid
        if self.__uid is None:
            self.__uid = ["uid"]
        self.__constraints = constraints
        self.__constraints_unique = unique_constraints
        self.__nonunique_index = nonunique_index

    def __hash__(self):
        return hash((self.neo4j_label))

    def _neo4j_repr(
        self,
        var="n",
    ):
        """Create a string of the form (var:LABEL {uid: 'uid', ...})"""
        uids = {i: self._Neo4jElement__properties.get(i) for i in self._Node__uid}
        uids_as_str = ", ".join(
            f"{k}: '{v}'" if isinstance(v, str) else f"{k}: {v}"
            for k, v in uids.items()
        )
        return f"({var}: {self._Neo4jElement__neo4j_label} {{{uids_as_str}}})"

    def _neo4j_repr_params(
        self,
        var="n",
        map_key="required_props",
    ):
        """Return a string of the form (var:LABEL {uid: $uid, ...})
        Ensure that the required properties are included in the string
        """
        req_props = self._Neo4jElement__required_properties_dict.keys()
        req_props = [f"{i}: {map_key}.{i}" for i in req_props]
        req_props = ", ".join(req_props)
        return f"({var}: {self._Neo4jElement__neo4j_label} {{{req_props}}})"

    def add_to_neo4j(
        self,
    ):
        """Add a single node to Neo4j"""

        merge_str = self._neo4j_repr_params(var="n", map_key="required_props")
        paramsdict = {
            "optional_props": self._Neo4jElement__optional_properties_dict,
            "required_props": self._Neo4jElement__required_properties_dict,
        }
        with GraphDriver() as db:
            results = db.run(
                f"""
                WITH $paramsdict as paramsdict
                WITH paramsdict.optional_props as optional_props, paramsdict.required_props as required_props
                MERGE {merge_str}
                ON CREATE SET n += optional_props
                """,
                paramsdict=paramsdict,
            ).value()

    @staticmethod
    def add_multiple_to_neo4j(list_of_nodes, batch_size=1000):
        """Add multiple nodes to Neo4j at a time in transactions of batch_size nodes at a time"""
        if not all(
            isinstance(sub, type(list_of_nodes[0])) for sub in list_of_nodes[1:]
        ):
            raise ValueError(
                f"All nodes in list_of_nodes must be of the same type as the first element, found (not in order): {set([sub.__class__.__name__ for sub in list_of_nodes])}"
            )
        single = list_of_nodes[0]
        merge_str = single._neo4j_repr_params(var="n", map_key="required_props")
        for paramsdictlist in batched(
            (
                {
                    "optional_props": i._Neo4jElement__optional_properties_dict,
                    "required_props": i._Neo4jElement__required_properties_dict,
                }
                for i in list_of_nodes
            ),
            batch_size,
        ):
            with GraphDriver() as db:
                results = db.run(
                    f"""
                    WITH $paramsdictlist as paramsdictlist
                    UNWIND paramsdictlist as paramsdict
                    WITH paramsdict.optional_props as optional_props, paramsdict.required_props as required_props
                    MERGE {merge_str}
                    ON CREATE SET n += optional_props
                    """,
                    paramsdictlist=paramsdictlist,
                ).value()


class Relationship(Neo4jElement):
    def __init__(self, start=None, end=None, **kwargs):
        super().__init__(**kwargs)
        self.start = start
        self.end = end
        if not self.start and not self.end:
            if self._Neo4jElement__header is not None:
                self.start = self._extract_start_from_admin_header
                self.end = self._extract_end_from_admin_header

    def _start_field(self):
        for i in self._Neo4jElement__header:
            if "START_ID" in i:
                return i

    def _end_field(self):
        for i in self._Neo4jElement__header:
            if "END_ID" in i:
                return i

    @property
    def _extract_start_from_admin_header(self):
        try:
            return re.findall(r"\((.*?)\)", self._start_field())[0]
        except Exception as e:
            raise e

    @property
    def _extract_end_from_admin_header(self):
        try:
            return re.findall(r"\((.*?)\)", self._end_field())[0]
        except Exception as e:
            raise e



    @property
    def _cypher_string(self):
        if self.start and self.end:
            return f"(:{self.start._Neo4jElement__neo4j_label})-[:{self._Neo4jElement__neo4j_label}]->(:{self.end._Neo4jElement__neo4j_label})"

    def add_to_neo4j(self):
        match_start = self.start._neo4j_repr_params(var="m1", map_key="required_props1")
        match_end = self.end._neo4j_repr_params(var="m2", map_key="required_props2")
        self.start.add_to_neo4j()
        self.end.add_to_neo4j()
        the_json = {
            "required_props1": self.start._Neo4jElement__required_properties_dict,
            "required_props2": self.end._Neo4jElement__required_properties_dict,
            "properties": self._Neo4jElement__properties,
        }
        with GraphDriver() as db:
            results = db.run(
                f"""
                WITH $properties as properties
                WITH properties.required_props1 as required_props1, properties.required_props2 as required_props2, properties.properties as props
                MATCH {match_start}
                MATCH {match_end}
                MERGE (m1)-[r:{self._Neo4jElement__neo4j_label}]->(m2)
                ON CREATE SET r += props
                """,
                properties=the_json,
            ).value()


