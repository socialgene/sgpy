import re
from abc import ABC
from itertools import batched
from typing import List

from socialgene.neo4j.neo4j import GraphDriver
from socialgene.utils.logging import log


class Neo4jElement(ABC):
    neo4j_label = None  # list for nodes because they can have multiple labels, string for relationships because they can only have one
    description = None
    property_specification = {}
    header = None
    header_filename = None
    target_subdirectory = None
    target_extension = None
    required_properties = None

    def __init__(
        self,
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
        """
        # properties aren't required for relationships
        if self.property_specification is None:
            self.property_specification = {}
        # if not provided then all properties are required

        if self.required_properties is None:
            self.required_properties = list(self.property_specification.keys())
        if self.required_properties is None:
            self.required_properties = []
        if properties is None:
            self.properties = {}
        if properties is not None:
            self.properties = properties
            self._clean_properties()

    @property
    def __required_properties_dict(self):
        return {i: self.properties.get(i) for i in self.required_properties}

    @property
    def __optional_properties_dict(self):
        return {
            i: self.properties.get(i)
            for i in self.properties
            if i not in self.required_properties
        }

    def __hash__(self):
        return hash((self.neo4j_label, frozenset(self.properties.items())))

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            if self.neo4j_label == other.neo4j_label:
                if self.properties == other.properties:
                    return True
        return False

    def _check_required_properties(
        self,
    ):
        missing = [i for i in self.required_properties if i not in self.properties]
        if missing:
            raise ValueError(
                f"{self.__class__} was missing required properties: {missing}"
            )

    def _check_no_extra_properties(self):
        extra = [i for i in self.properties if i not in self.property_specification]
        if extra:
            raise ValueError(
                f"{self.__class__} contained unsupported extra properties: {extra}"
            )

    def _check_property_types(self):
        bad = {
            k: f"expected {self.property_specification.get(k).__name__} but got {type(v).__name__}"
            for k, v in self.properties.items()
            if not isinstance(v, self.property_specification.get(k))
        }
        if bad:
            raise ValueError(
                f"{self.__class__} contained properties with incorrect types: {bad}"
            )
        self.properties = {
            k: v
            for k, v in self.properties.items()
            if isinstance(v, self.property_specification.get(k))
        }

    def _clean_properties(
        self,
    ):
        # only keep properties that are in the property specification and remove nulls
        self.properties = {
            k: v
            for k, v in self.properties.items()
            if k in self.property_specification and v is not None
        }
        self._check_required_properties()
        self._check_no_extra_properties()
        self._check_property_types()

    def fill_from_dict(self, d: dict):
        """Fill the properties of the object from a dictionary"""
        for k, v in self.property_specification.items():
            try:
                if k not in d:
                    log.debug(
                        f"Failed to fill add property {k} from dictionary: not found for {self.__class__}"
                    )
                    continue
                if d.get(k) is None:
                    log.debug(
                        f"Failed to fill add property {k} from dictionary: None for {self.__class__}"
                    )
                    continue

                self.properties[k] = v(d.get(k))
            except Exception as e:
                log.debug(
                    f"Failed to fill add property {k} from dictionary: {e} for {self.__class__}"
                )


class Node(Neo4jElement):
    """Represents a single Node"""

    nonunique_index = None

    def __init__(
        self,
        uid: List[str] = None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        if self.property_specification is None:
            raise ValueError(
                f"property_specification must be defined for class {self.__class__.__name__}"
            )

    def __hash__(self):
        return hash((frozenset(self.neo4j_label), frozenset(self.properties.items())))

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            if self.neo4j_label == other.neo4j_label:
                if self.properties == other.properties:
                    return True
        return False

    def __str__(self) -> str:
        return f'(:{":".join(self.neo4j_label)})'

    def _neo4j_repr(
        self,
        var="n",
    ):
        """Create a string of the form (var:LABEL {uid: 'uid', ...})"""
        uids = {i: self.properties.get(i) for i in self._Node__uid}
        uids_as_str = ", ".join(
            f"{k}: '{v}'" if isinstance(v, str) else f"{k}: {v}"
            for k, v in uids.items()
        )
        return f"({var}: {self.neo4j_label} {{{uids_as_str}}})"

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
        return f"({var}: {self.neo4j_label[0]} {{{req_props}}})"

    def add_constraints_to_neo4j(self):
        """Add constraints to Neo4j"""
        if not hasattr(self, "constraints_unique") or not self.constraints_unique:
            log.info(f"No unique constraints to add for {self.__str__()} nodes")
            return
        # for now at least only use the first label
        label = self.neo4j_label[0]
        with GraphDriver() as db:
            res = db.run(
                f"""
                CREATE CONSTRAINT python_added_index_{label} IF NOT EXISTS
                FOR (n:{label})
                REQUIRE ({", ".join([f"n.{i}" for i in self.constraints_unique])}) IS UNIQUE
                """
            ).consume()
            try:
                # The code is checking if the variable `res` has a property called `counters` and if
                # it does, it will execute the code block inside the `if` statement.
                if res.counters:
                    log.info(
                        f"Created {res.counters.constraints_added} constraints for {self.__str__()} nodes"
                    )
            except Exception:
                pass

    def add_nonunique_index_to_neo4j(self):
        """Add non-unique indices to Neo4j"""
        if not hasattr(self, "nonunique_index") or not self.nonunique_index:
            log.info(f"No non-unique indices to add for {self.__str__()} nodes")
            return
        label = self.neo4j_label[0]
        with GraphDriver() as db:
            res = db.run(
                f"""
                CREATE INDEX python_added_non_unique_index_{label} IF NOT EXISTS
                FOR (n:{label})
                ON ({", ".join([f"n.{i}" for i in self.nonunique_index])})
                """
            ).consume()
            try:
                if res.counters:
                    log.info(
                        f"Created {res.counters.indexes_added} non-unique indices for {self.__str__()} nodes"
                    )
            except Exception:
                pass

    def add_to_neo4j(self, create=False):
        """Add a single node to Neo4j"""
        var = "n"
        merge_str = self._neo4j_repr_params(var=var, map_key="required_props")
        paramsdict = {
            "optional_props": self._Neo4jElement__optional_properties_dict,
            "required_props": self._Neo4jElement__required_properties_dict,
        }
        if len(self.neo4j_label) > 1:
            # handle adding extra labels if they exist
            add_label_str = ":".join(self.neo4j_label[1:])
            add_label_str = f"SET {var}:{add_label_str}"
        else:
            add_label_str = ""
        if create:
            with GraphDriver() as db:
                res = db.run(
                    f"""
                    WITH $paramsdict as paramsdict
                    WITH paramsdict.optional_props as optional_props, paramsdict.required_props as required_props
                    CREATE {merge_str}
                    SET n += optional_props
                    {add_label_str}
                    """,
                    paramsdict=paramsdict,
                ).consume()
                try:
                    log.info(
                        f"Created {res.counters.nodes_created} {self.__str__()} nodes, set {res.counters.properties_set} properties"
                    )
                except Exception:
                    pass
        else:
            if len(self.neo4j_label) > 1:
                add_label_str = ":".join(self.neo4j_label[1:])
                add_label_str = f"ON CREATE SET {var}:{add_label_str}"
            else:
                add_label_str = ""
            with GraphDriver() as db:
                res = db.run(
                    f"""
                    WITH $paramsdict as paramsdict
                    WITH paramsdict.optional_props as optional_props, paramsdict.required_props as required_props
                    MERGE {merge_str}
                    ON CREATE SET n += optional_props
                    {add_label_str}
                    """,
                    paramsdict=paramsdict,
                ).consume()
                try:
                    log.info(
                        f"Created {res.counters.nodes_created} {self.__str__()} nodes, set {res.counters.properties_set} properties"
                    )
                except Exception:
                    pass

    @staticmethod
    def add_multiple_to_neo4j(list_of_nodes, batch_size=1000, create=False):
        """Add multiple nodes to Neo4j at a time in transactions of batch_size nodes at a time"""
        if not all(
            isinstance(sub, type(list_of_nodes[0])) for sub in list_of_nodes[1:]
        ):
            raise ValueError(
                f"All nodes in list_of_nodes must be of the same type as the first element, found (not in order): {set([sub.__class__.__name__ for sub in list_of_nodes])}"
            )
        single = list_of_nodes[0]
        count_res = {"nodes_created": 0, "properties_set": 0}
        var = "n"
        merge_str = single._neo4j_repr_params(var=var, map_key="required_props")
        if len(single.neo4j_label) > 1:
            # handle adding extra labels if they exist
            add_label_str = ":".join(single.neo4j_label[1:])
            if create:
                add_label_str = f"SET {var}:{add_label_str}"
            else:
                add_label_str = f"ON CREATE SET {var}:{add_label_str}"
        else:
            add_label_str = ""
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
            if create:
                with GraphDriver() as db:
                    res = db.run(
                        f"""
                        WITH $paramsdictlist as paramsdictlist
                        UNWIND paramsdictlist as paramsdict
                        WITH paramsdict.optional_props as optional_props, paramsdict.required_props as required_props
                        CREATE {merge_str}
                        SET n += optional_props
                        {add_label_str}
                        """,
                        paramsdictlist=paramsdictlist,
                    ).consume()
                try:
                    count_res["nodes_created"] += res.counters.nodes_created
                    count_res["properties_set"] += res.counters.properties_set
                except Exception:
                    pass
            else:
                with GraphDriver() as db:
                    res = db.run(
                        f"""
                        WITH $paramsdictlist as paramsdictlist
                        UNWIND paramsdictlist as paramsdict
                        WITH paramsdict.optional_props as optional_props, paramsdict.required_props as required_props
                        MERGE {merge_str}
                        ON CREATE SET n += optional_props
                        {add_label_str}
                        """,
                        paramsdictlist=paramsdictlist,
                    ).consume()
                try:
                    count_res["nodes_created"] += res.counters.nodes_created
                    count_res["properties_set"] += res.counters.properties_set
                except Exception:
                    pass
        log.info(
            f"Created {count_res['nodes_created']} {single.__str__()} nodes, set {count_res['properties_set']} properties"
        )


class Relationship(Neo4jElement):
    def __init__(self, start=None, end=None, **kwargs):
        super().__init__(**kwargs)
        self.start = start
        self.end = end
        if not isinstance(self.start, self.start_class):
            raise ValueError(
                f"Relationship class '{self.__class__.__name__}' start node must be of type {self.start_class.__name__}, found {self.start.__class__.__name__}"
            )
        if not isinstance(self.end, self.end_class):
            raise ValueError(
                f"Relationship class '{self.__class__.__name__}' end node must be of type {self.end_class.__name__}, found {self.end.__class__.__name__}"
            )

    def __str__(self):
        return f"{self.start.__str__()}-[:{self.neo4j_label}]->{self.end.__str__()}"

    def __hash__(self):
        return hash(
            (
                self.neo4j_label,
                frozenset(self.properties.items()),
                self.start,
                self.end,
            )
        )

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            if self.neo4j_label == other.neo4j_label:
                if self.properties == other.properties:
                    if self.start == other.start:
                        if self.end == other.end:
                            return True
        return False

    @property
    def neo4j_label_str(self):
        return str(self.neo4j_label)

    def _start_field(self):
        for i in self.header:
            if "START_ID" in i:
                return i

    def _end_field(self):
        for i in self.header:
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
        return f"(:{self.start_class.neo4j_label})-[:{self.neo4j_label}]->(:{self.end_class.neo4j_label})"

    def add_to_neo4j(self, create=False):
        match_start = self.start._neo4j_repr_params(var="m1", map_key="required_props1")
        match_end = self.end._neo4j_repr_params(var="m2", map_key="required_props2")
        # self.start.add_to_neo4j()
        # self.end.add_to_neo4j()
        the_json = {
            "required_props1": self.start._Neo4jElement__required_properties_dict,
            "required_props2": self.end._Neo4jElement__required_properties_dict,
            "properties": self.properties,
        }
        if create:
            with GraphDriver() as db:
                x = db.run(
                    f"""
                    WITH $properties as properties
                    WITH properties.required_props1 as required_props1, properties.required_props2 as required_props2, properties.properties as props
                    MATCH {match_start}
                    MATCH {match_end}
                    CREATE (m1)-[r:{self.neo4j_label}]->(m2)
                    CALL apoc.create.relationship(m1, '{self.neo4j_label}', props, m2)
                    YIELD rel
                    RETURN count(rel) as relationships_created
                    """,
                    properties=the_json,
                ).value()
        else:
            with GraphDriver() as db:
                x = db.run(
                    f"""
                    WITH $properties as properties
                    WITH properties.required_props1 as required_props1, properties.required_props2 as required_props2, properties.properties as props
                    MATCH {match_start}
                    MATCH {match_end}
                    CALL apoc.merge.relationship(m1, '{self.neo4j_label}', props, {{}}, m2, {{}})
                    YIELD rel
                    RETURN count(rel) as relationships_created
                    """,
                    properties=the_json,
                ).value()
        log.info(f"{x[0]} relationships created {self.__str__()}")

    @staticmethod
    def add_multiple_to_neo4j(list_of_rels, batch_size=1000, create=False):
        """Add multiple relationships to Neo4j at a time in transactions of batch_size nodes at a time"""
        if not all(isinstance(sub, type(list_of_rels[0])) for sub in list_of_rels[1:]):
            raise ValueError(
                f"All relationships in list_of_rels must be of the same type as the first element, found (not in order): {set([sub.__class__.__name__ for sub in list_of_rels])}"
            )
        single = list_of_rels[0]
        match_start = single.start._neo4j_repr_params(
            var="m1", map_key="required_props1"
        )
        match_end = single.end._neo4j_repr_params(var="m2", map_key="required_props2")
        n_rels_created = 0
        for paramsdictlist in batched(
            (
                {
                    "required_props1": i.start._Neo4jElement__required_properties_dict,
                    "required_props2": i.end._Neo4jElement__required_properties_dict,
                    "properties": i.properties,
                }
                for i in list_of_rels
            ),
            batch_size,
        ):
            if create:
                with GraphDriver() as db:
                    x = db.run(
                        f"""
                        WITH $paramsdictlist as paramsdictlist
                        UNWIND paramsdictlist as properties
                        WITH properties.required_props1 as required_props1, properties.required_props2 as required_props2, properties.properties as props
                        MATCH {match_start}
                        MATCH {match_end}
                        CALL apoc.create.relationship(m1, '{single.neo4j_label}', props, m2 )
                        YIELD rel
                        RETURN count(rel) as relationships_created
                        """,
                        paramsdictlist=paramsdictlist,
                    ).value()
                n_rels_created += x[0]
            else:
                with GraphDriver() as db:
                    x = db.run(
                        f"""
                        WITH $paramsdictlist as paramsdictlist
                        UNWIND paramsdictlist as properties
                        WITH properties.required_props1 as required_props1, properties.required_props2 as required_props2, properties.properties as props
                        MATCH {match_start}
                        MATCH {match_end}
                        CALL apoc.merge.relationship(m1, '{single.neo4j_label}', props, {{}}, m2, {{}})
                        YIELD rel
                        RETURN count(rel) as relationships_created
                        """,
                        paramsdictlist=paramsdictlist,
                    ).value()
                n_rels_created += x[0]
        log.info(f"{n_rels_created} relationships created {single.__str__()}")
