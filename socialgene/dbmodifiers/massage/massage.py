import json
from pathlib import Path
import tempfile
import requests
import tarfile
from socialgene.dbmodifiers.mibig.nrps import NRPS
from socialgene.neo4j.neo4j import GraphDriver
from socialgene.utils.logging import log
from socialgene.config import env_vars
from socialgene.dbmodifiers.mibig.compound import Mibig_Compound
from rich import inspect
from neo4j import GraphDatabase
from socialgene.config import env_vars
import logging
from rich.progress import Progress, SpinnerColumn
from rich.console import Console
import time

console = Console()


def _add_to_neo4j(statement, **kwargs):
    try:
        summary = (
            GraphDriver().driver.execute_query(
                statement,
                **kwargs,
                database_="neo4j",
            )
        ).summary
        if summary.metadata.get("stats"):
            log.info(
                f"{summary.metadata.get('stats').get('properties-set')} properties modified"
            )
        else:
            log.info("No properties modified")
    except Exception as e:
        log.warning(e.message)


def _run_transaction_function(statement, **kwargs):
    try:
        with GraphDriver() as db:
            a = db.run(
                statement,
                **kwargs,
                database_="neo4j",
            ).values()
        return a

    except Exception as e:
        log.error(e)


class Massage:
    @staticmethod
    def add_taxonomic_name_to_assembly(rank):
        _add_to_neo4j(
            f"""
            MATCH (a1:assembly)-[:IS_TAXON|TAXON_PARENT*1..]->(t1:taxid {{rank:"{rank}"}})
            SET a1.{rank} = t1.name;
            """
        )

    @staticmethod
    def add_protein_descriptions():
        res = _run_transaction_function(
            """
                // set protein description as a random input protein's description (except "hypothetical")
                MATCH (p1:protein)
                CALL {
                    WITH p1
                    MATCH ()-[e1:ENCODES]->(p1:protein)
                    WHERE e1.description IS NOT NULL AND NOT e1.description CONTAINS "hypothetical"
                    RETURN e1 LIMIT 1
                } IN TRANSACTIONS OF 1000 ROWS
                SET p1.description = e1.description
                RETURN count(p1)
            """
        )
        log.info(f"Modified: {res} properties")

    @staticmethod
    def add_antismash_regions():
        _run_transaction_function(
            """
                CALL apoc.load.json("import/antismash_results.jsonl")
                    YIELD value AS jsonl
                CALL {
                    with jsonl
                    UNWIND jsonl.records as record
                    UNWIND keys(record) as nucleotide_id
                    WITH record, nucleotide_id, range(0, size(record[nucleotide_id])-1) as cluster_index
                    UNWIND cluster_index as cluster_index2
                    WITH record, record[nucleotide_id][cluster_index2] as cluster, nucleotide_id,  cluster_index2
                    MATCH (n1:nucleotide {external_id:nucleotide_id})-[c1:ENCODES]->()
                    WHERE cluster.start < (c1.start + 1 ) < cluster.end
                    SET c1.antismash_products = record[nucleotide_id][cluster_index2].products
                    SET c1.antismash_region = cluster_index2
                } IN TRANSACTIONS OF 1000 ROWS;
            """
        )

    @staticmethod
    def culture_collections_as_nodes_rels():
        _add_to_neo4j(
            """
                MATCH (a1:assembly) where a1.culture_collection is not null
                WITH a1, split(a1.culture_collection,":")[0] as cc_agency
                MERGE (cc:culture_collection {uid:cc_agency})
                MERGE (a1)-[:FOUND_IN]->(cc);
            """
        )

    @staticmethod
    def set_mibig_bgc():
        _run_transaction_function(
            """
                MATCH (a1:assembly)<-[:ASSEMBLES_TO]-(n1:nucleotide)
                WHERE a1.uid STARTS WITH "BGC00" AND n1.external_id STARTS WITH "BGC00"
                SET a1:mibig_bgc;
            """
        )
