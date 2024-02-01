import builtins
from textwrap import wrap

from rich import print
from rich.console import Console, ConsoleOptions, RenderResult
from rich.table import Table
from socialgene.addons.chebi import Chebi
from socialgene.addons.chembl import Chembl
from socialgene.addons.classyfire import ClassyFire
from socialgene.addons.npatlas import NPAtlasEntry
from socialgene.addons.npmrd import Npmrd

from socialgene.base.molbio import LocusAssemblyMetadata
from socialgene.config import env_vars
from socialgene.neo4j.neo4j_element import Neo4jElement, Node
from socialgene.utils.lists_to_markdown import markdown_table_from_list

# use rich to print
builtins.print = print


class NodesMixin:
    """Represents multiple Nodes and is where all addon Nodes are defined"""

    def __init__(
        self,
    ):
        self.nodes = {}

        self.add_node(
            neo4j_label="npatlas",
            description="Represents an NPAtlas entry.",
            header=sorted(NPAtlasEntry.__slots__),
        )
        self.add_node(
            neo4j_label="chebi",
            description="Represents a Chebi entry.",
            header=sorted(Chebi.__slots__),
        )
        self.add_node(
            neo4j_label="chembl",
            description="Represents a Chembl entry.",
            header=sorted(Chembl.__slots__),
        )
        self.add_node(
            neo4j_label="classyfire",
            description="Represents a ClassyFire entry.",
            header=sorted(ClassyFire.__slots__),
        )
        self.add_node(
            neo4j_label="npmrd",
            description="Represents a Chebi entry.",
            header=sorted(Npmrd.__slots__),
        )

    def add_node(self, neo4j_label, **kwargs):
        if neo4j_label in self.nodes:
            raise ValueError(f"Node with label {neo4j_label} already exists")
        self.nodes[neo4j_label] = Node(neo4j_label=neo4j_label, **kwargs)

    def __rich_console__(
        self, console: Console, options: ConsoleOptions
    ) -> RenderResult:  # pragma: no cover
        table = Table(title="Nodes", show_lines=True)
        table.add_column("Label", justify="left", style="cyan", no_wrap=True, ratio=1)
        table.add_column(
            "Description",
            justify="left",
            style="cyan",
            no_wrap=False,
            ratio=4,
            max_width=50,
        )
        table.add_column("Properties", style="magenta", ratio=1)
        # sort by label which is the key
        for i in (self.nodes[i] for i in sorted(self.nodes.keys())):
            table.add_row(
                i._Neo4jElement__neo4j_label,
                i._Neo4jElement__description,
                "\n".join(wrap(", ".join(i._Neo4jElement__header))),
            )
        yield table

    def _markdown_table(self):
        cols = [
            (
                "Label",
                "Description",
                "NF results subdirectory",
                "Neo4j header file",
                "header",
            )
        ]
        rows = [
            (
                i.neo4j_label,
                i.description,
                i.target_subdirectory,
                i.header_filename,
                i.header,
            )
            for i in (self.nodes[i] for i in sorted(self.nodes.keys()))
        ]
        cols.extend(rows)
        print(
            markdown_table_from_list(
                cols,
                align="left",
            )
        )


def print_info():  # pragma: no cover
    print(NodesMixin())


def print_markdown():  # pragma: no cover
    NodesMixin()._markdown_table()


def printer():  # pragma: no cover
    import argparse

    parser = argparse.ArgumentParser(description="Print node info")
    parser.add_argument(
        "--markdown",
        help="",
        default=False,
        required=False,
        action=argparse.BooleanOptionalAction,
    )
    args = parser.parse_args()
    if args.markdown:
        print_markdown()
    else:
        print_info()


if __name__ == "__main__":
    print_info()
