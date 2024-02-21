from rdkit import DataStructs

from socialgene.addons.chemistry.nr import ChemicalCompoundNode, ChemicalSimilarity
from socialgene.base.chem import ChemicalCompound
from socialgene.neo4j.neo4j import GraphDriver
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    TextColumn,
    TimeElapsedColumn,
)
import multiprocessing
import pandas as pd
from itertools import batched
cmpd_label = ChemicalCompoundNode.neo4j_label


with GraphDriver() as db:
    res = db.run(
        f"""
        MATCH (c1:{cmpd_label})
        RETURN c1.inchi as chem
        """,
    ).value()


def make_node(x):
    res = {}
    for i in x:
        temp = ChemicalCompound(i)
        res[i] = (temp, temp.node)
    return res



num_processes = multiprocessing.cpu_count() - 2
group_size = len(res) // num_processes

b = batched(res, group_size)
combined_dict = {}
with multiprocessing.Pool(processes=num_processes) as pool:
    for res in pool.imap_unordered(make_node, b):
        combined_dict.update(res)




id_list = list(combined_dict.keys())

big_vec = [i[0].morgan for i in combined_dict.values()]
cutoff = 0.5




def calculate_similarity(tupgroup):
    chemsim_set = set()
    tup = tupgroup[2]
    for iter in range(tupgroup[0], tupgroup[1]):
        for i, x in enumerate(DataStructs.BulkTanimotoSimilarity(tup[iter], tup)):
            if x >= 0.5:
                if tup[iter] != i:
                    score = int(x * 100)
                    temp = (iter, i, score)
                    chemsim_set.add(temp)
    return chemsim_set


# Split x into groups
x=[(i,big_vec) for i in range(len(id_list))]

group_size = len(id_list) // num_processes
groups = [(i,i+group_size, big_vec) for i in range(0, len(id_list), group_size)]
groups[-1] = (groups[-1][0], len(id_list) - 1, big_vec)
chemsim = set()

with multiprocessing.Pool(processes=num_processes) as pool:
    for i in pool.imap_unordered(calculate_similarity, groups):
        chemsim.update(i)


chemsim2 = set()

for i in chemsim:
    # sorted prevents identical forward + reverse relationships
    chemsim2.add(tuple(sorted(i[0:2]) + [i[2]]))

chemsim=chemsim2
del chemsim2

# get nodes faster using pd indexing
df0 = pd.DataFrame(((k, v[0], v[1]) for k, v in combined_dict.items()), columns = ["id", "compound", "node"])
df = pd.DataFrame(chemsim, columns=["start", "end", "similarity"])
df = df.set_index('start')
df = df.join(df0['node'])
df.columns = ["end", "similarity", "start_node"]
df = df.reset_index(level=0)
df = df.set_index('end')
df = df.join(df0['node'])
df.columns = ["start", "similarity", "start_node", "end_node"]
df = df.reset_index(level=0)
df['relationship'] = df.apply(lambda x: ChemicalSimilarity(start=x.start_node, end=x.end_node, properties={"similarity": x.similarity}), axis=1)


zz={i:z[id_list[i]].node for i,x in enumerate(id_list)}

df['relationship'] = df.apply(lambda x: ChemicalSimilarity(start=x.start_node.node, end=x.end_node.node, properties={"similarity": x.similarity}), axis=1)


a=[ChemicalSimilarity(start=z[id_list[i[0]]].node,end=z[id_list[i[1]]].node,properties={"similarity": i[2]},) for i in chemsim]
ChemicalSimilarity(start=z[id_list[i[0]]].node,end=z[id_list[i[1]]].node,properties={"similarity": i[2]},)




