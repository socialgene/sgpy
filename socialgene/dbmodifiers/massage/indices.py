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


def _add_to_neo4j(statement, **kwargs):
    try:
        (
            GraphDriver().driver.execute_query(
                statement,
                **kwargs,
                database_="neo4j",
            )
        )
    except Exception as e:
        # If fails because index exists, don't fail, just report it. Fail any other errors.
        if e.code == "Neo.ClientError.Schema.EquivalentSchemaRuleAlreadyExists":
            log.warning(e.message)
        else:
            raise e


class Indices:
    def mess(func):
        def wrapper():
            log.info("May appear to hang while index/constraint is populating")
            func()

        return wrapper

    @staticmethod
    def constraint_status():
        with GraphDriver() as db:
            print(db.run("SHOW CONSTRAINTS").to_df())

    @staticmethod
    def index_status():
        with GraphDriver() as db:
            print(db.run("SHOW INDEXES").to_df())

    @mess
    @staticmethod
    def assembly_uid():
        _add_to_neo4j(
            """
            CREATE CONSTRAINT assembly_uid IF NOT EXISTS
            FOR (n:assembly)
            REQUIRE (n.uid) IS UNIQUE;
            """
        )

    @mess
    @staticmethod
    def protein_uid():
        _add_to_neo4j(
            """
            CREATE CONSTRAINT protein_uid IF NOT EXISTS
            FOR (n:protein)
            REQUIRE (n.uid) IS UNIQUE;
            """
        )

    @mess
    @staticmethod
    def hmm_uid():
        _add_to_neo4j(
            """
            CREATE CONSTRAINT hmm_uid IF NOT EXISTS
            FOR (n:hmm)
            REQUIRE (n.uid) IS UNIQUE;
            """
        )

    @mess
    @staticmethod
    def nucleotide_uid():
        _add_to_neo4j(
            """
            CREATE CONSTRAINT nucleotide_uid IF NOT EXISTS
            FOR (n:nucleotide)
            REQUIRE (n.uid) IS UNIQUE;
            """
        )

    @mess
    @staticmethod
    def nucleotide_external_id():
        _add_to_neo4j(
            """
            CREATE INDEX nuc_extern_id
            FOR (n:nucleotide)
            ON (n.external_id);
            """
        )
