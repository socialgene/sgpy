"""https://mibig.secondarymetabolites.org/"""
import re

from socialgene.addons.base import ExternalBaseClass


class Mibig(ExternalBaseClass):
    __slots__ = ["uid"]

    def __init__(self, uid) -> None:
        self.uid = self.attempt_to_fix_uid(uid)

    @staticmethod
    def attempt_to_fix_uid(uid):
        # npatlas has some mibig ids prefixed, some not, try to add prefix here
        # https://github.com/NPAtlas/npatlas_website/issues/13
        if uid.startswith("BGC"):
            temp = uid
        else:
            temp = "BGC" + uid.zfill(7)
        # Validate uid
        if re.match("BGC[0-9]{7}", temp):
            return temp
        else:
            raise ValueError(f"Unexpected mibig ID {uid}")

    def add_node_to_neo4j(self):
        # empty because this should connect to mibig assembly nodes
        pass
