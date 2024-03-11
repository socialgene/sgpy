# flake8: noqa F401
import builtins
from textwrap import wrap

from rich.table import Table

import socialgene.addons.chebi.nr
import socialgene.addons.chembl.nr
import socialgene.addons.chemistry.nr
import socialgene.addons.classyfire.nr
import socialgene.addons.gnps_library.nr
import socialgene.addons.gnps_networking.nr
import socialgene.addons.mibig.nr
import socialgene.addons.npatlas.nr
import socialgene.addons.npclassifier.nr
import socialgene.addons.npmrd.nr
import socialgene.addons.publication.nr
import socialgene.nextflow.nodes
import socialgene.nextflow.relationships
from socialgene.config import env_vars
from socialgene.neo4j.neo4j_element import Node, Relationship
from socialgene.utils.lists_to_markdown import markdown_table_from_list
from socialgene.utils.logging import CONSOLE, log

builtins.print = print


def recursive_get(x, levels=10):
    if levels == 0:
        return set(x)
    else:
        subclasses = [i.__subclasses__() for i in x]
        flattened_subclasses = [
            subclass for sublist in subclasses for subclass in sublist
        ]
        return set(x) | recursive_get(flattened_subclasses, levels - 1)


class GraphSchema:
    ALL_NODES = recursive_get(Node.__subclasses__())
    ALL_RELATIONSHIPS = recursive_get(Relationship.__subclasses__())
    NEXTFLOW_NODES = [
        x for x in ALL_NODES if x.__module__.startswith("socialgene.nextflow")
    ]
    NEXTFLOW_RELATIONSHIPS = [
        x for x in ALL_RELATIONSHIPS if x.__module__.startswith("socialgene.nextflow")
    ]
    ADDON_NODES = [x for x in ALL_NODES if x.__module__.startswith("socialgene.addons")]
    ADDON_RELATIONSHIPS = [
        x for x in ALL_RELATIONSHIPS if x.__module__.startswith("socialgene.addons")
    ]

    def __init__(self):
        pass

    def _nodes_table(self):
        table = Table(title="Nodes", show_lines=True)
        table.add_column("Label", justify="left", style="cyan", no_wrap=True, ratio=1)
        table.add_column("From", justify="left", style="cyan", no_wrap=True, ratio=1)
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
        for i in sorted(list(self.ALL_NODES), key=lambda x: x.neo4j_label[0]):
            table.add_row(
                ":".join(i.neo4j_label),
                i.__module__,
                i.description,
                "\n".join(wrap(", ".join(i.property_specification))),
            )
        yield table

    def _relationships_table(self):
        table = Table(title="Relationships", show_lines=True)
        table.add_column("Label", justify="left", style="cyan", no_wrap=True, ratio=1)
        table.add_column(
            "Defined in", justify="left", style="cyan", no_wrap=True, ratio=1
        )
        table.add_column("Relationship", style="magenta", ratio=1)
        table.add_column("NF results subdirectory", style="magenta", ratio=1)
        table.add_column("Neo4j header file", style="magenta", ratio=1)
        for i in sorted(list(self.ALL_RELATIONSHIPS), key=lambda x: x.neo4j_label[0]):
            table.add_row(
                i.neo4j_label,
                i.__module__,
                i(i.start_class(), i.end_class()).__str__(),
                i.target_subdirectory,
                i.header_filename,
            )
        yield table

    def _markdown_table_nodes(nodelist):
        cols = [
            (
                "Label",
                "Description",
                "NF results subdirectory",
                "Neo4j header file",
                "Unique on",
                "properties",
            )
        ]
        rows = [
            (
                ":".join(i.neo4j_label),
                i.description,
                i.target_subdirectory,
                i.header_filename,
                ", ".join(i.constraints_unique),
                [k for k in i.property_specification.keys()],
            )
            for i in sorted(list(nodelist), key=lambda x: x.neo4j_label[0])
        ]
        cols.extend(rows)
        print(
            markdown_table_from_list(
                cols,
                align="left",
            )
        )

    def _markdown_table_rels(rellist):
        cols = [
            (
                "Label",
                "Relationship",
                "NF results subdirectory",
                "Neo4j header file",
            )
        ]
        rows = [
            (
                i.neo4j_label,
                i(i.start_class(), i.end_class()).__str__(),
                i.target_subdirectory,
                i.header_filename,
            )
            for i in sorted(list(rellist), key=lambda x: x.neo4j_label[0])
        ]
        cols.extend(rows)

        print(
            markdown_table_from_list(
                cols,
                align="left",
            )
        )


def main():  # pragma: no cover
    import argparse

    parser = argparse.ArgumentParser(description="Print node info")
    parser.add_argument(
        "--markdown",
        help="",
        default=False,
        required=False,
        action=argparse.BooleanOptionalAction,
    )
    parser.add_argument(
        "--nodes",
        help="",
        default=False,
        required=False,
        action=argparse.BooleanOptionalAction,
    )
    parser.add_argument(
        "--rels",
        help="",
        default=False,
        required=False,
        action=argparse.BooleanOptionalAction,
    )
    args = parser.parse_args()
    if args.markdown:
        if args.nodes:
            GraphSchema._markdown_table_nodes(GraphSchema.ALL_NODES)
        if args.rels:
            GraphSchema._markdown_table_rels(GraphSchema.ALL_RELATIONSHIPS)
    else:
        before_width = CONSOLE.width
        CONSOLE.width = 300
        if args.nodes:
            CONSOLE.print(GraphSchema()._nodes_table().__next__())
            CONSOLE.width = before_width
        if args.rels:
            CONSOLE.print(GraphSchema()._relationships_table().__next__())
        CONSOLE.width = before_width


if __name__ == "__main__":
    main()
