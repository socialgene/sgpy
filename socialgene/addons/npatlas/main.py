import argparse
from collections import defaultdict
from socialgene.addons.classyfire.parse import ClassyFire
from socialgene.addons.npatlas.parse import NPAtlas, NPAtlasEntry
import json
from rich.progress import Progress
from rich.progress import Progress, BarColumn, TimeElapsedColumn,TextColumn
import concurrent.futures

a = NPAtlas(atlas_json_path="/home/chase/Downloads/NPAtlas_download.json")


nodes = set()
rel_nodes = defaultdict(set)
rel_rels = defaultdict(set)
cc = 0

def process_entry(entry):
    z = NPAtlasEntry(entry)
    z.parse()
    nodes.add(z.node)
    for k, v in z.rel_rels.items():
        rel_rels[k].update(v)
    for k, v in z.rel_nodes.items():
        rel_nodes[k].update(v)



with Progress(
        TextColumn("{task.completed}"),
        "[progress.description]{task.description}",
        BarColumn(),
        "[progress.percentage]{task.percentage:>3.0f}%",
        TimeElapsedColumn(),
        ) as pg:
    task = pg.add_task("[cyan]Processing NPAtlas entries...", total=33372)
    with open(a.path) as f:
        entries = json.load(f)  # Limit the number of entries to process
        with concurrent.futures.ThreadPoolExecutor() as executor:
            for entry in concurrent.futures.as_completed([executor.submit(process_entry, entry) for entry in entries]):
                pg.update(task, advance=1)


list(nodes)[0].add_multiple_to_neo4j(list(nodes), create=True)


f=ClassyFire()
f.add_classyfire_nodes()
f.add_classyfire_isa_relationships()
f.add_classyfire_synonym_relationships()



for k,v in rel_nodes.items():
    print(f"adding {k} nodes")
    list(v)[0].add_multiple_to_neo4j(list(v), create=False)

for k,v in rel_rels.items():
    print(f"adding {k} relationships")
    list(v)[0].add_multiple_to_neo4j(list(v), create=False)



