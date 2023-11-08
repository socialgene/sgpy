import argparse
import csv
import re
from copy import deepcopy
from pathlib import Path

import requests

from socialgene.utils.logging import log

parser = argparse.ArgumentParser(
    description="Download and create edgelist of all Go terms"
)

parser.add_argument(
    "--outdir",
    help="Output directory filepath",
    required=True,
)

OBO_URL = "http://purl.obolibrary.org/obo/go/go-basic.obo"

goterm_pattern = re.compile("GO:[0-9]{7}", flags=0)


def extract_goterm_int_as_str(id):
    if goterm := goterm_pattern.search(id):
        return str(goterm.group().removeprefix("GO:"))


class GoNode:
    __slots__ = ["uid", "name", "namespace", "definition"]

    def __init__(self) -> None:
        self.uid = None
        self.name = None
        self.namespace = None
        self.definition = None

    def add_id(self, x):
        if goterm := extract_goterm_int_as_str(x):
            self.uid = goterm
        else:
            raise ValueError(f"'{x}' isn't a goterm")

    def add_name(self, x):
        self.name = x.removeprefix("name:").strip()

    def add_namespace(self, x):
        self.namespace = x.removeprefix("namespace:").strip()

    def add_definition(self, x):
        self.definition = x.removeprefix("definition:").strip()

    def assign(self, line):
        if line.startswith("id:"):
            self.add_id(line)
        elif line.startswith("name:"):
            self.add_name(line)
        elif line.startswith("namespace:"):
            self.add_namespace(line)
        elif line.startswith("definition:"):
            self.add_definition(line)

    def output(self):
        return (self.uid, self.name, self.namespace, self.definition)


class GoRelationship:
    __slots__ = [
        "start",
        "end",
        "type",
    ]

    def __init__(self) -> None:
        self.start = None
        self.end = None
        self.type = None

    def add_end(self, id):
        if goterm := extract_goterm_int_as_str(id):
            self.end = goterm
        else:
            print(f"'{id}' isn't a goterm")

    def add_type(self, line):
        key = line.split(":", 1)[0]
        if key == "relationship":
            # 'relationship: part_of GO:0000785 ! chromatin'
            # should return: "part of"
            self.type = line.split(":", 1)[1].strip().split(" ", 1)[0]
        elif key == "is_a":
            self.type = key
        self.type = self.type.upper()
        self.end = extract_goterm_int_as_str(line)

    def assign(self, start, line):
        self.start = start
        if line.startswith("is_a") or line.startswith("relationship"):
            self.add_type(line)

    def output(self):
        if all((self.start, self.end, self.type)):
            return (self.start, self.end, self.type)


class AltRel(GoRelationship):
    def __init__(
        self,
    ):
        super().__init__()

    def assign(self, main_id, line):
        if line.startswith("alt_id"):
            self.start = main_id
            self.end = extract_goterm_int_as_str(line)
            self.type = "ALTERNATE"


def write(
    input: list[tuple[str, ...]],
    outpath: str,
) -> None:
    with open(outpath, "w") as h:
        writer = csv.writer(
            h, quoting=csv.QUOTE_MINIMAL, delimiter="\t", dialect="unix"
        )
        writer.writerows(input)


def download_obo(url=OBO_URL):
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        for line in r.iter_lines():
            yield line.decode(r.encoding).strip()


def read_from_file(filepath):
    with open(filepath, "r") as h:
        for i in h:
            yield i.strip()


def parse(outdir, filepath=None, url=None):
    nodes = []
    relationships = []
    if url:
        _open = download_obo(url=url)
    elif filepath:
        _open = read_from_file(filepath=filepath)
    else:
        raise ValueError("Must input url or filepath")
    term_open = False
    for line in _open:
        if line == "[Term]":
            term_open = True
            nodeobj = GoNode()
            relobj = GoRelationship()
        if not term_open:
            continue
        nodeobj.assign(line=line)
        # some goterms have "alt_ids" these are sometimes referenced from other sources but only exist
        # in the obo file as alt ids (no primary entry)
        if line.startswith("alt_id"):
            # deepcopy so name and namespace will carry if present
            alt_node = deepcopy(nodeobj)
            alt_node.add_id(line)
            alt_rel = AltRel()
            alt_rel.assign(main_id=nodeobj.uid, line=line)
            nodes.append(alt_node.output())
            relationships.append(alt_rel.output())
            del alt_node, alt_rel
        relobj.assign(start=nodeobj.uid, line=line)
        if line == "":
            term_open = False
            nodes.append(nodeobj.output())
            # only output full relationhip tuples
            if out := relobj.output():
                relationships.append(out)
    log.info("Writing goterms nodes")
    write(nodes, Path(outdir, "goterms"))
    log.info("Writing goterms relationships")
    write(relationships, Path(outdir, "goterm_edgelist"))
    log.info("Finished reading/writing goterms")


def main():
    args = parser.parse_args()
    parse(url=OBO_URL, outdir=args.outdir)


if __name__ == "__main__":
    main()
