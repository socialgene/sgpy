from typing import List

from socialgene.neo4j.neo4j import Neo4jQuery

neo4j_object = Neo4jQuery()


def search_protein_hash(query_hash: List[str]):
    """Search the database for proteins with exact hash matches

    Args:
        query_hash (list): list of hash strings

    Raises:
        ValueError: incorrect input

    Returns:
        dict: dictionary of whether the protein hash is in the database, e.g. {query1hash: True, query2hash: False}
    """
    if isinstance(query_hash, str):
        query_hash = [query_hash]
    elif not isinstance(query_hash, list) and all(
        [isinstance(i, str) for i in query_hash]
    ):
        raise ValueError("query_hash must be a list of strings")
    # results will be e.g. [{'result': ['XcXHqKUXe3T2zDq3bhJxuw_0rJ2zwreU']}]
    results = neo4j_object.query_neo4j(
        cypher_name="search_protein_hash", param=query_hash
    )
    # Create a dictionary of all input values, with False as the default
    query_hash_dict = {i: False for i in query_hash}
    # Change returned ids to True
    for i in results[0]["result"]:
        query_hash_dict[i] = True
    return query_hash_dict
