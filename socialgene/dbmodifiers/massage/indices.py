from socialgene.neo4j.neo4j import GraphDriver
from socialgene.utils.logging import log


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


def mess(func):
    def wrapper():
        log.info("May appear to hang while index/constraint is populating")
        func()

    return wrapper


def constraint_status():
    with GraphDriver() as db:
        print(db.run("SHOW CONSTRAINTS").to_df())


def index_status():
    with GraphDriver() as db:
        print(db.run("SHOW INDEXES").to_df())


@mess
def assembly_uid():
    _add_to_neo4j(
        """
            CREATE CONSTRAINT assembly_uid IF NOT EXISTS
            FOR (n:assembly)
            REQUIRE (n.uid) IS UNIQUE;
            """
    )


@mess
def protein_uid():
    _add_to_neo4j(
        """
            CREATE CONSTRAINT protein_uid IF NOT EXISTS
            FOR (n:protein)
            REQUIRE (n.uid) IS UNIQUE;
            """
    )


@mess
def hmm_uid():
    _add_to_neo4j(
        """
            CREATE CONSTRAINT hmm_uid IF NOT EXISTS
            FOR (n:hmm)
            REQUIRE (n.uid) IS UNIQUE;
            """
    )


@mess
def nucleotide_uid():
    _add_to_neo4j(
        """
            CREATE CONSTRAINT nucleotide_uid IF NOT EXISTS
            FOR (n:nucleotide)
            REQUIRE (n.uid) IS UNIQUE;
            """
    )


@mess
def nucleotide_external_id():
    _add_to_neo4j(
        """
            CREATE INDEX nuc_extern_id
            FOR (n:nucleotide)
            ON (n.external_id);
            """
    )


@mess
def taxonomy_uid():
    _add_to_neo4j(
        """
            CREATE INDEX taxonomy_uid
            FOR (n:taxid)
            ON (n.uid);
            """
    )


@mess
def gnps_library_uid():
    _add_to_neo4j(
        """
            CREATE INDEX gnps_library_uid
            FOR (n:gnps_library_spectrum)
            ON (n.uid);
            """
    )


@mess
def npatlas_uid():
    _add_to_neo4j(
        """
            CREATE INDEX npatlas_uid
            FOR (n:npatlas)
            ON (n.uid);
            """
    )
