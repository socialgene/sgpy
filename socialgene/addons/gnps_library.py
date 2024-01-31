import re
from typing import List

from socialgene.addons.base import ExternalBaseClass


class GnpsLibrarySpectrum(ExternalBaseClass):
    __slots__ = ["uid"]

    def __init__(self, uid) -> None:
        self.uid = uid

    @staticmethod
    def _extract_CCMSLIB(x) -> List:
        return re.findall("CCMSLIB[0-9]{11}", x)

    def add_node_to_neo4j(self):
        self._add_to_neo4j(
            """
            MERGE (:gnps {uid: $uid})
            """,
            uid=self.uid,
        )
