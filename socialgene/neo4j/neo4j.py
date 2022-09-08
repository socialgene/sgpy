# python dependencies
import os
import importlib.resources
import atexit
from typing import Any

# external dependencies
from neo4j import GraphDatabase

# internal dependencies
from socialgene.utils.logging import log
from socialgene.config import env_vars

# TODO: Not sure whether the driver should be persistent (with active check)
# or a new connection each time


# Read in queries from cypher file, as dictionary
def import_queries():
    cypher_dictionary = {}
    with importlib.resources.open_text("socialgene.neo4j", "queries.cypher") as f:
        for line in f:
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


# class Neo4j:
#     __slots__ = ["driver", "address"]
#     neo4jVersion = os.getenv("NEO4J_VERSION", "4")
#     username = os.getenv("NEO4J_USER", "neo4j")
#     password = os.getenv("NEO4J_PASSWORD", "test")
#     neo4jVersion = os.getenv("NEO4J_VERSION", "4")
#     database = os.getenv("NEO4J_DATABASE", "neo4j")
#     port = os.getenv("PORT", 8080)

#     def __init__(self):
#         self.driver = None
#         # for non-django, run this before starting python:
#         # export NEO4J_URI="neo4j://localhost:7687"
#         self.address = env_vars["NEO4J_URI"]

#     def __enter__(self):
#         self.driver = GraphDatabase.driver(
#             self.address, auth=(self.username, self.password)
#         )
#         try:
#             self.driver.verify_connectivity()
#         except Exception as e:
#             log.error(e)
#             self.driver.verify_connectivity()
#         return self.driver.session()

#     def __exit__(self, exc_type, exc_value, exc_traceback):
#         self.driver.close()


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
                        os.getenv("NEO4J_USER", "neo4j"),
                        os.getenv("NEO4J_PASSWORD", "test"),
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
    def query_neo4j(cypher_name: str, param: Any):
        """Run a provided cypher query on {param}

        Args:
            cypher_name (str): Neo4j Cypher query
            param (Any): paramter passed to Neo4j query, type depends on query

        Returns:
            Any: type depends on query
        """
        # See queries.py
        resolved_query_dict = import_queries()
        # grab the the neo4j connection
        with GraphDriver() as db:
            # make the query against the db
            results = db.read_transaction(
                lambda tx: tx.run(
                    resolved_query_dict[cypher_name]["query"], param=param
                ).data()
            )
        return results
