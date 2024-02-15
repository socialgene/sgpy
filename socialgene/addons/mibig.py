"""https://mibig.secondarymetabolites.org/"""

import re

from socialgene.neo4j.neo4j_element import Node
from socialgene.nextflow.nodes import ASSEMBLY


class Mibig(ASSEMBLY):
    # TODO: multiple labels
    neo4j_label = "assembly"
    description = "Represents a single Mibig entry"
    property_specification = {
        "uid": str,
    }
    required_properties = ["uid"]

    def __init__(self, properties) -> None:
        super().__init__()
        self.properties["uid"] = self.attempt_to_fix_uid(properties["uid"])

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

