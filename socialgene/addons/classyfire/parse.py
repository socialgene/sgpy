"""http://classyfire.wishartlab.com"""

import re
from io import BytesIO
from zipfile import ZipFile

import requests

from socialgene.addons.chebi.nr import ChebiNode
from socialgene.addons.classyfire.nr import (
    ClassyFireIsA,
    ClassyFireNode,
    ClassyFireSynonym,
)
from socialgene.utils.logging import log

OBO_URL = (
    "http://classyfire.wishartlab.com/system/downloads/1_0/chemont/ChemOnt_2_1.obo.zip"
)


class ClassyFireEntry:
    __slots__ = ["uid", "name", "definition", "chebi", "is_a"]

    def __init__(self) -> None:
        self.uid = None
        self.name = None
        self.definition = None
        self.chebi = set()
        self.is_a = None

    def assign(self, line):
        if line.startswith("id:"):
            self._add_id(line)
        elif line.startswith("name:"):
            self._add_name(line)
        elif line.startswith("def:"):
            self._add_definition(line)
        elif line.startswith("synonym:"):
            self._add_chebilist(line)
        elif line.startswith("is_a:"):
            self._add_is_a(line)

    def _add_id(self, x):
        self.uid = x.removeprefix("id: CHEMONTID:").strip()
        self.uid = int(self.uid)

    def _add_name(self, x):
        self.name = x.removeprefix("name:").strip()

    def _add_definition(self, x):
        self.definition = x.split('"')[1].strip()

    def _add_chebilist(self, x):
        matched = re.search("CHEBI:[0-9]+", x)
        if matched:
            self.chebi.add(matched.group().removeprefix("CHEBI:").strip())

    def _add_is_a(self, x):
        matched = re.search("CHEMONTID:[0-9]+", x)
        if matched:
            self.is_a = matched.group().removeprefix("CHEMONTID:").strip()
            self.is_a = int(self.is_a)


class ClassyFire:
    def __init__(self) -> None:
        self.entries = self.download()

    @staticmethod
    def download(url=OBO_URL):
        response = requests.get(url, stream=True)
        obs = []
        log.info(f"Downloading and parsing data from:\n\t{OBO_URL}")
        with ZipFile(BytesIO(response.content)).open("ChemOnt_2_1.obo") as h:
            term_open = False
            for line in h:
                line = line.decode("utf-8")
                if line.startswith("[Term]"):
                    term_open = True
                    obj = ClassyFireEntry()
                if not term_open:
                    continue
                obj.assign(line)
                if line == "\n":
                    term_open = False
                    obs.append(obj)
        return obs

    def add_classyfire_nodes(self, create=False):
        b = [
            ClassyFireNode(
                properties={"uid": i.uid, "name": i.name, "definition": i.definition}
            )
            for i in self.entries
        ]
        ClassyFireNode.add_multiple_to_neo4j(b, create=create)

    def add_classyfire_isa_relationships(self, create=False):
        b = [
            ClassyFireIsA(
                start=ClassyFireNode(properties={"uid": i.uid}),
                end=ClassyFireNode(properties={"uid": i.is_a}),
                create=create,
            )
            for i in self.entries
            if i.is_a
        ]
        ClassyFireIsA.add_multiple_to_neo4j(b, create=create)

    def add_classyfire_synonym_relationships(self, create=False):
        zz = []
        chebi_ids = set()
        for i in self.entries:
            for ii in i.chebi:
                chebi_ids.add(ii)
                zz.append(
                    ClassyFireSynonym(
                        start=ClassyFireNode(properties={"uid": i.uid}),
                        end=ChebiNode(properties={"uid": int(ii)}),
                        create=create,
                    )
                )
        ChebiNode.add_multiple_to_neo4j(
            [ChebiNode(properties={"uid": int(i)}) for i in chebi_ids], create=create
        )
        ClassyFireSynonym.add_multiple_to_neo4j(zz, create=create)
