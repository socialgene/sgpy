import argparse
import concurrent.futures
import json
import tempfile
from collections import defaultdict
from pathlib import Path

from rich.progress import BarColumn, Progress, TextColumn, TimeElapsedColumn

from socialgene.addons.npatlas.parse import NPAtlasEntry, _download_npatlas
from socialgene.utils.logging import log


def create_arg_parser():
    """ "Creates and returns the ArgumentParser object."""
    parser = argparse.ArgumentParser(
        description="Integrate NPAtlas into a SocialGene Neo4j Database"
    )
    parser.add_argument(
        "--input",
        help="Local path to NPAtlas JSON file. If None, downloads from NPAtlas website.",
        default=None,
        required=False,
    )
    return parser


def main():
    parser = create_arg_parser()
    args = parser.parse_args()
    input = args.input
    if input is None:
        with tempfile.NamedTemporaryFile() as temp:
            Path(temp.name).exists()
            _download_npatlas(outpath=temp.name)
            with open(temp.name) as f:
                entries = json.load(f)
    else:
        jpath = Path(input)
        if jpath.exists():
            with open(jpath) as f:
                entries = json.load(f)
        else:
            raise FileNotFoundError(f"Path {input} does not exist")

    nodes = set()
    rels = defaultdict(set)

    def process_entry(entry):
        z = NPAtlasEntry(entry)
        z.parse()
        return z

    with Progress(
        TextColumn("{task.completed}"),
        "[progress.description]{task.description}",
        BarColumn(),
        "[progress.percentage]{task.percentage:>3.0f}%",
        TimeElapsedColumn(),
    ) as pg:
        task = pg.add_task("[cyan]Processing NPAtlas entries...", total=33372)
        # This was multithreaded whcih worked well on the first iteraion of this code when it was primarily reading
        # but after chem data calculation was added, it may be just as slow as single
        with concurrent.futures.ThreadPoolExecutor() as executor:
            for entry in concurrent.futures.as_completed(
                [executor.submit(process_entry, entry) for entry in entries]
            ):
                z = entry.result()
                nodes.add(z.node)
                for k, v in z.get_links().items():
                    rels[k].update(v)
                pg.update(task, advance=1)
    log.info("Creating/Merging npatlas nodes in neo4")
    list(nodes)[0].add_multiple_to_neo4j(list(nodes), create=False)
    relnodes = defaultdict(set)
    for k, v in rels.items():
        for i in v:
            if i.start.__class__.__name__ != "NPAtlasNode":
                relnodes[i.start.__class__].add(i.start)
            if i.end.__class__.__name__ != "NPAtlasNode":
                relnodes[i.end.__class__].add(i.end)

    log.info("Creating/Merging nodes linked to npatlas entries in neo4j")
    # make sure nodes that npatlas entries are linked to exist in the database
    for k, v in relnodes.items():
        list(v)[0].add_multiple_to_neo4j(list(v), create=False)
    log.info("Linking npatlas entries and related nodes in neo4j")
    for k, v in rels.items():
        k.add_multiple_to_neo4j(list(v), create=False)


if __name__ == "__main__":
    main()
