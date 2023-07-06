import re
from socialgene.neo4j.neo4j import GraphDriver


class ClassyFire:
    def __init__(self) -> None:
        self.uid = None
        self.name = None
        self.definition = None
        self.chebi = set()
        self.is_a = None

    def assign(self, line):
        if line.startswith("id:"):
            self.add_id(line)
        elif line.startswith("name:"):
            self.add_name(line)
        elif line.startswith("def:"):
            self.add_definition(line)
        elif line.startswith("synonym:"):
            self.add_chebilist(line)
        elif line.startswith("is_a:"):
            self.add_is_a(line)

    def add_id(self, x):
        self.uid = x.removeprefix("id: CHEMONTID:").strip()

    def add_name(self, x):
        self.name = x.removeprefix("name:").strip()

    def add_definition(self, x):
        self.definition = x.split('"')[1].strip()

    def add_chebilist(self, x):
        matched = re.search("CHEBI:[0-9]+", x)
        if matched:
            self.chebi.add(matched.group().removeprefix("CHEBI:").strip())

    def add_is_a(self, x):
        matched = re.search("CHEMONTID:[0-9]+", x)
        if matched:
            self.is_a = matched.group().removeprefix("CHEMONTID:").strip()

    def _create_chemont_nodes(self):
        _ = GraphDriver().driver.execute_query(
            """
                MERGE (:chemont { uid: $uid, name: $name, definition: $definition } )
                """,
            uid=self.uid,
            name=self.name,
            definition=self.definition,
            database_="neo4j",
        )

    def create_chebi_nodes(self):
        _ = GraphDriver().driver.execute_query(
            """
            WITH $chebi as chebi_uids
            UNWIND chebi_uids as chebi_uid
            MERGE ( :chebi { uid: chebi_uid} )
            """,
            chebi=list(self.chebi),
            database_="neo4j",
        )

    def connect_chemont_chebi_nodes(self):
        _ = GraphDriver().driver.execute_query(
            """
            WITH $chebi as chebi_uids
            UNWIND chebi_uids as chebi_uid
            MATCH (chem:chemont {uid:$uid})
            MATCH (chebi:chebi {uid: chebi_uid})
            MERGE (chem)-[:SYNONYM]->(chebi)

            """,
            uid=self.uid,
            chebi=list(self.chebi),
            database_="neo4j",
        )

    def connect_chemont_is_a_nodes(self):
        _ = GraphDriver().driver.execute_query(
            """
            MATCH (chem1:chemont {uid:$uid})
            MATCH (chem2:chemont {uid:$is_a_uid})
            MERGE (chem1)-[:IS_A]->(chem2)

            """,
            uid=self.uid,
            is_a_uid=self.is_a,
            database_="neo4j",
        )
