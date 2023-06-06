import json
from pathlib import Path
from rich import inspect
from rdkit import Chem
from rdkit import DataStructs
from socialgene.neo4j.neo4j import GraphDriver  # grab the the neo4j connection
import logging
from rdkit.Chem import AllChem
from rich.progress import Progress
import pandas as pd
import itertools

json_dir = "/home/chase/Documents/data/mibig/3_1/mibig_json_3.1"

pathlist = Path(json_dir).glob("BGC*.json")


collect_dict = {}

for json_path in pathlist:
    try:
        with open(json_path, "r") as f:
            data = json.load(f)
        if "cluster" in data:
            if "compounds" in data["cluster"]:
                # only get the first one
                collect_dict[data["cluster"]["mibig_accession"]] = data["cluster"][
                    "compounds"
                ][0]["chem_struct"]
    except:
        pass


fpgen = AllChem.GetRDKitFPGenerator()
m = {k: Chem.MolFromSmiles(v) for k, v in collect_dict.items()}

fps = {k: fpgen.GetFingerprint(v) for k, v in m.items()}

z = []
with Progress(transient=True) as pg:
    task = pg.add_task("Progress...", total=1690041)
    for i in itertools.combinations(fps, 2):
        z.append((i[0], i[1], DataStructs.TanimotoSimilarity(fps[i[0]], fps[i[1]])))
    pg.update(task, advance=1)


df = pd.DataFrame(z)


logging.getLogger("neo4j").setLevel(logging.WARNING)
logging.getLogger().setLevel(logging.INFO)

query = """
WITH $params as inputs
UNWIND inputs AS input
MATCH (a1:assembly)
WHERE a1.id = input[0]
SET a1.SMILES = input[1]
"""

with GraphDriver() as db:
    results = db.run(query, params=[[k, v] for k, v in collect_dict.items()])


query = """
WITH $params as inputs
UNWIND inputs AS input
WITH input
WHERE input[2] > 0.9
MATCH (a2:assembly {id: input[0]})
MATCH (a1:assembly {id: input[1]})
CREATE (a1)-[:TANIMOTO_SIMILARITY {score: input[2]}]->(a2)

"""

with GraphDriver() as db:
    results = db.run(query, params=z)


with GraphDriver() as db:
    with Progress(transient=True) as pg:
        task = pg.add_task("Progress...", total=len([i for i in z if i[2] > 0.5]))
        for i in [i for i in z if i[2] > 0.5]:
            results = db.run(
                """
        WITH $params as input
        MATCH (a2:assembly {id: input[0]})
        MATCH (a1:assembly {id: input[1]})
        CREATE (a1)-[:TANIMOTO_SIMILARITY {score: input[2]}]->(a2)
        """,
                params=i,
            )
        pg.update(task, advance=1)


fpgen = AllChem.GetRDKitFPGenerator()

m = {k: Chem.MolFromSmiles(v) for k, v in collect_dict.items()}

fps = {k: fpgen.GetFingerprint(v) for k, v in m.items()}

z = []
with Progress(transient=True) as pg:
    task = pg.add_task("Progress...", total=1690041)
    for i in itertools.combinations(fps, 2):
        z.append((i[0], i[1], DataStructs.TanimotoSimilarity(fps[i[0]], fps[i[1]])))
    pg.update(task, advance=1)


COCONUT_DB.smi
