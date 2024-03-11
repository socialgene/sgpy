import argparse

from rich.progress import Progress, SpinnerColumn, TextColumn

from socialgene.neo4j.schema.graph_schema import GraphSchema
from socialgene.utils.logging import CONSOLE, log

parser = argparse.ArgumentParser(
    description="Add indices to a SocialGene Neo4j Database"
)
parser.add_argument(
    "--all",
    help="",
    default=False,
    required=False,
    action="store_true",
)
parser.add_argument(
    "--labels",
    help="Node labels to try and add indices to",
    default=False,
    required=False,
    nargs="+",
)


def main():
    args = parser.parse_args()
    as_dict = {}
    log.info(
        "Adding indices to Neo4j. This can take some time, depending on how many of each label exists in the database."
    )
    for i in GraphSchema.ALL_NODES:
        if len(i.neo4j_label) == 1:
            as_dict[i.neo4j_label[0]] = i
        else:
            as_dict[i.neo4j_label[-1]] = i
    if not any([i for i in args.__dict__.values()]):
        parser.print_help()
        print(f"Available labels: {sorted(list(as_dict.keys()))}")
        return
    if args.all:
        with Progress(
            SpinnerColumn(spinner_name="runner"),
            TextColumn(text_format="Adding index for... {task.fields[lab]}"),
            console=CONSOLE,
            transient=True,
        ) as progress:
            task = progress.add_task("Adding index for...", lab="")
            for i in GraphSchema.ALL_NODES:
                progress.update(task, lab=i.neo4j_label)
                try:
                    i().add_constraints_to_neo4j()
                    i().add_nonunique_index_to_neo4j()
                except Exception as e:
                    log.warning(f"Failed to add indices to {i.neo4j_label}: {e}")

    elif args.labels:
        with Progress(
            SpinnerColumn(spinner_name="runner"),
            TextColumn(text_format="Adding index for... {task.fields[lab]}"),
            console=CONSOLE,
            transient=True,
        ) as progress:
            task = progress.add_task("Adding index for...", lab="")
            for i in args.labels:
                if i in as_dict:
                    progress.update(task, lab=i)
                    try:
                        as_dict[i]().add_constraints_to_neo4j()
                        as_dict[i]().add_nonunique_index_to_neo4j()
                    except Exception as e:
                        log.warning(f"Failed to add indices to {i}: {e}")
                else:
                    log.warning(f"Label {i} not found in GraphSchema.ALL_NODES")
                    log.warning(f"Available labels: {as_dict.keys()}")


if __name__ == "__main__":
    main()
