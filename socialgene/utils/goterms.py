import argparse
import csv
from pathlib import Path
import requests
import re

from socialgene.utils.logging import log

parser = argparse.ArgumentParser(
    description="Download and create edgelist of all Go terms"
)

parser.add_argument(
    "--outdir",
    metavar="filepath",
    help="Output directory filepath",
    required=True,
)

OBO_URL = "http://purl.obolibrary.org/obo/go/go-basic.obo"

goterm_pattern = re.compile("GO:[0-9]{7}", flags=0)


def extract_goterm_int_as_str(id):
    if goterm := goterm_pattern.search(id):
        return str(goterm.group().removeprefix("GO:"))


class GoNode:
    __slots__ = ["id", "name", "namespace", "definition"]

    def __init__(self) -> None:
        self.id = None
        self.name = None
        self.namespace = None
        self.definition = None

    def add_id(self, x):
        if goterm := extract_goterm_int_as_str(x):
            self.id = goterm
        else:
            print(f"'{x}' isn't a goterm")
            raise

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
        return (self.id, self.name, self.namespace, self.definition)


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
        self.type = f"GO_{self.type.upper()}"
        self.end = extract_goterm_int_as_str(line)

    def assign(self, start, line):
        self.start = start
        if line.startswith("is_a") or line.startswith("relationship"):
            self.add_type(line)

    def output(self):
        return (self.start, self.end, self.type)


def write(
    input: list[tuple[str, ...]],
    outpath: str,
) -> None:
    with open(outpath, "w") as h:
        writer = csv.writer(
            h, quoting=csv.QUOTE_MINIMAL, delimiter="\t", dialect="unix"
        )
        writer.writerows(input)


def main():
    args = parser.parse_args()
    nodes = []
    relationships = []
    with requests.get(OBO_URL, stream=True) as r:
        r.raise_for_status()
        term_open = False
        for line in r.iter_lines():
            line = line.decode(r.encoding)
            if line == "[Term]":
                term_open = True
                nodeobj = GoNode()
                relobj = GoRelationship()
            if not term_open:
                # skip lines before first term
                continue
            nodeobj.assign(line=line)
            relobj.assign(start=nodeobj.id, line=line)
            if line == "":
                term_open = False
                nodes.append(nodeobj.output())
                relationships.append(relobj.output())

    log.info("Writing goterms nodes")
    write(nodes, Path(args.outdir, "goterms"))
    log.info("Writing goterms relationships")
    write(relationships, Path(args.outdir, "goterm_edgelist"))
    log.info("Finished reading/writing goterms")


if __name__ == "__main__":
    main()
