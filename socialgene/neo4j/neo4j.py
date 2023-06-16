from typing import Any

import atexit
import importlib.resources

from neo4j import GraphDatabase

from socialgene.config import env_vars
from socialgene.utils.logging import log

# TODO: Not sure whether the driver should be persistent (with active check)
# or a new connection each time


# Read in queries from cypher file, as dictionary
def import_queries():
    cypher_dictionary = {}
    with open(
        importlib.resources.files("socialgene").joinpath("neo4j", "queries.cypher"), "r"
    ) as f:
        for line in f:
            pass
            if line.startswith("// Name:"):
                query_name = line.replace("// Name:", "").strip()
                cypher_dictionary[query_name] = {
                    "query": "",
                    "description": "",
                    "param": "",
                }
            elif line.startswith("// Description:"):
                temp = line.replace("// Description:", "").strip()
                cypher_dictionary[query_name]["description"] = temp
                del temp
            elif line.startswith("// Param:"):
                temp = line.replace("// Param:", "").strip()
                cypher_dictionary[query_name]["Param"] = temp
                del temp
            else:
                cypher_dictionary[query_name][
                    "query"
                ] = f"{cypher_dictionary[query_name]['query']}\n{line.strip()}"
    return cypher_dictionary


class GraphDriver(object):
    __instance = None
    __driver = None

    def __new__(cls):
        if cls.__instance is None:
            inst = cls.__instance = object.__new__(cls)
            try:
                inst.__driver = GraphDatabase.driver(
                    env_vars["NEO4J_URI"],
                    auth=(
                        env_vars["NEO4J_USER"],
                        env_vars["NEO4J_PASSWORD"],
                    ),
                )
                log.info(f"Connected to Neo4j database at {env_vars['NEO4J_URI']}")
            except Exception as e:
                log.exception("Error connecting to Neo4j database")
                log.exception(e)
        atexit.register(cls.shutdown)
        return cls.__instance

    def __init__(self) -> None:
        self.session = None

    def check_connection(self):
        self.__driver.verify_connectivity()

    def __enter__(self):
        self.session = self.__driver.session()
        return self.session

    def __exit__(self, *args, **kwargs):
        self.session.close()

    @classmethod
    def shutdown(cls):
        if cls.__instance:
            cls.__instance.__driver.close()
            cls.__instance = None


class Neo4jQuery:
    def __init__(self):
        pass

    @staticmethod
    def print_query(query_name):  # pragma: no cover
        log.info(import_queries()[query_name]["query"])

    @staticmethod
    def print_queries():  # pragma: no cover
        resolved_query_dict = import_queries()
        for k in resolved_query_dict.keys():
            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            log.info(f"Query name: {k}")
            log.info(f'Description: \t{resolved_query_dict[k]["description"]}')
            log.info(f'Parameter(s): \t{resolved_query_dict[k]["param"]}')
            log.info(f'Query: {resolved_query_dict[k]["query"]}')
            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

    @staticmethod
    def query_neo4j(
        cypher_name: str = None,
        cypher: str = None,
        param: Any = None,
        rettype="data",
        *args,
        **kwargs,
    ):
        """Run a provided cypher query on {param}

        Args:
            cypher_name (str, optional): Neo4j Cypher query. Defaults to None.
            cypher (str, optional): Single string of Cypher script.  Defaults to None.
            param (Any, optional): paramter passed to Neo4j query, type depends on query. Defaults to None.
            rettype (str, optional): output function. For available methods see (https://neo4j.com/docs/api/python-driver/current/api.html#result). Defaults to "data"
            args/kwargs (Any): pas additional argument(s) to the rettype method
        Returns:
            Any: type depends on query
        """
        # See queries.py
        resolved_query_dict = import_queries()
        # grab the the neo4j connection
        if cypher_name:
            query = resolved_query_dict[cypher_name]["query"]
        elif cypher:
            query = cypher
        with GraphDriver() as db:
            # make the query against the db
            results = db.run(query, param=param)
            return getattr(results, rettype)(*args, **kwargs)
