from textwrap import wrap
import socialgene.nextflow.nodes
import socialgene.nextflow.relationships

from socialgene.neo4j.neo4j_element import Node, Relationship
from socialgene.utils.lists_to_markdown import markdown_table_from_list
from rich.console import Console, ConsoleOptions, RenderResult
from rich.table import Table
import builtins

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

    def __rich_console__(
        self, console: Console, options: ConsoleOptions
    ) -> RenderResult:  # pragma: no cover
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
                i._Neo4jElement__neo4j_label,
                i._Neo4jElement__description,
                "\n".join(wrap(", ".join(i._Neo4jElement__header))),
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


def print_info():  # pragma: no cover
    print(GraphSchema().__rich_console__())


def print_markdown():  # pragma: no cover
    _markdown_table_nodes(1)


builtins.print = print


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
