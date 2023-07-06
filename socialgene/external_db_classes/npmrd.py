from socialgene.external_db_classes.base_class import ExternalBaseClass


class Npmrd:
    __slots__ = ["uid"]

    def __init__(self, uid) -> None:
        self.uid = uid

    def add_node_to_neo4j(self):
        self._add_to_neo4j(
            """
            MERGE (:npmrd {uid: $uid})
            """,
            uid=self.uid,
        )
