from textwrap import wrap
import socialgene.nextflow.nodes
import socialgene.nextflow.relationships

from socialgene.neo4j.neo4j_element import Node, Relationship
from socialgene.utils.lists_to_markdown import markdown_table_from_list
from rich.console import Console, ConsoleOptions, RenderResult
from rich.table import Table
import builtins
from rich.console import Console

from socialgene.config import env_vars
from socialgene.utils.logging import log

console = Console()
builtins.print = print


class GraphSchema:
    ALL_NODES = Node.__subclasses__()
    ALL_RELATIONSHIPS = Relationship.__subclasses__()
    NEXTFLOW_NODES = [
        x for x in ALL_NODES if x.__module__.startswith("socialgene.nextflow")
    ]
    NEXTFLOW_RELATIONSHIPS = [
        x for x in ALL_RELATIONSHIPS if x.__module__.startswith("socialgene.nextflow")
    ]
    ADDON_NODES = [x for x in ALL_NODES if x.__module__.startswith("socialgene.addons")]

    def __init__(self):
        pass

    def _nodes_table(self):
        node_dict = {i.__name__: i for i in self.ALL_NODES}
        node_dict = dict(sorted(node_dict.items()))
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
        for i in node_dict.values():
            table.add_row(
                i()._Neo4jElement__neo4j_label,
                i()._Neo4jElement__description,
                "\n".join(wrap(", ".join(i()._Neo4jElement__header))),
            )
        yield table

    def _relationships_table(self):
        rel_dict = {i.__name__: i for i in self.ALL_RELATIONSHIPS}
        rel_dict = dict(sorted(rel_dict.items()))
        table = Table(title="Relationships", show_lines=True)
        table.add_column("Label", justify="left", style="cyan", no_wrap=True, ratio=1)
        table.add_column("Relationship", style="magenta", ratio=1)
        table.add_column("NF results subdirectory", style="magenta", ratio=1)
        table.add_column("Neo4j header file", style="magenta", ratio=1)
        for i in rel_dict.values():
            table.add_row(
                i()._Neo4jElement__neo4j_label,
                i()._cypher_string,
                i()._Neo4jElement__target_subdirectory,
                i()._Neo4jElement__header_filename,
            )
        yield table

    def _markdown_table_nodes(nodelist):
        node_dict = {i.__name__: i for i in nodelist}
        node_dict = dict(sorted(node_dict.items()))
        cols = [
            (
                "Label",
                "Description",
                "NF results subdirectory",
                "Neo4j header file",
                "properties",
            )
        ]
        rows = [
            (
                i()._Neo4jElement__neo4j_label,
                i()._Neo4jElement__description,
                i()._Neo4jElement__target_subdirectory,
                i()._Neo4jElement__header_filename,
                i()._Neo4jElement__properties,
            )
            for i in node_dict.values()
        ]
        cols.extend(rows)
        print(
            markdown_table_from_list(
                cols,
                align="left",
            )
        )

    def _markdown_table_rels(rellist):
        rel_dict = {i.__name__: i for i in rellist}
        rel_dict = dict(sorted(rel_dict.items()))
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
                i()._Neo4jElement__neo4j_label,
                i()._cypher_string,
                i()._Neo4jElement__target_subdirectory,
                i()._Neo4jElement__header_filename,
            )
            for i in rel_dict.values()
        ]
        cols.extend(rows)
        print(
            markdown_table_from_list(
                cols,
                align="left",
            )
        )


builtins.print = print


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
        if args.nodes:
            console.print(GraphSchema()._nodes_table().__next__())
        if args.rels:
            console.print(GraphSchema()._relationships_table().__next__())


if __name__ == "__main__":
    main()
