"""https://np-mrd.org/"""


from socialgene.neo4j.schema.relationship_mixins import ASSEMBLES_TO


class Npmrd:
    __slots__ = ["uid"]
    ...
    # __slots__ = ["uid"]

    # def __init__(self, uid) -> None:
    #     self.uid = uid

    # def add_node_to_neo4j(self):
    #     self._add_to_neo4j(
    #         """
    #         MERGE (:npmrd {uid: $uid})
    #         """,
    #         uid=self.uid,
    #     )
