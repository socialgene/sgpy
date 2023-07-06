from socialgene.external_db_classes.base_class import ExternalBaseClass


class Mibig(ExternalBaseClass):
    __slots__ = ["uid"]

    def __init__(self, uid) -> None:
        self.uid = uid

    def add_node_to_neo4j(self):
        # empty because this should connect to mibig assembly nodes
        pass
