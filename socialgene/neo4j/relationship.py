import builtins

from rich import print
from rich.console import Console, ConsoleOptions, RenderResult
from rich.table import Table

from socialgene.neo4j.neo4j_element import Neo4jElement, Relationship
from socialgene.utils.lists_to_markdown import markdown_table_from_list

# use rich to print
builtins.print = print


class RelationshipsMixin:
    def __init__(self, **kwargs):
        super(RelationshipsMixin, self).__init__()
        self.relationships = {}

    def add_relationship(self, neo4j_label, **kwargs):
        if neo4j_label in self.relationships:
            raise ValueError(f"Relationship {neo4j_label} already exists")
        self.relationships[neo4j_label] = Relationship(
            neo4j_label=neo4j_label, **kwargs
        )

    def __rich_console__(
        self, console: Console, options: ConsoleOptions
    ) -> RenderResult:
        table = Table(title="Relationships", show_lines=True)
        table.add_column("Label", justify="left", style="cyan", no_wrap=True, ratio=1)
        table.add_column("Relationship", style="magenta", ratio=1)
        table.add_column("NF results subdirectory", style="magenta", ratio=1)
        table.add_column("Neo4j header file", style="magenta", ratio=1)
        for i in (self.relationships[i] for i in sorted(self.relationships.keys())):
            table.add_row(
                i.neo4j_label,
                i.cypher_string,
                i.target_subdirectory,
                i.header_filename,
            )
        yield table

    def _markdown_table(self):
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
                i.cypher_string,
                i.target_subdirectory,
                i.header_filename,
            )
            for i in (self.relationships[i] for i in sorted(self.relationships.keys()))
        ]
        cols.extend(rows)
        print(
            markdown_table_from_list(
                cols,
                align="left",
            )
        )


def print_info():  # pragma: no cover
    print(RelationshipsMixin())


def print_markdown():  # pragma: no cover
    RelationshipsMixin()._markdown_table()


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
