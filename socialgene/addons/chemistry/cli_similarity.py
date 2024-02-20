import argparse
from socialgene.addons.chemistry.nr import ChemicalCompoundNode,ChemicalSimilarity
from socialgene.base.chem import ChemicalCompound
from rdkit import DataStructs

from socialgene.neo4j.neo4j import GraphDriver



cmpd_label = ChemicalCompoundNode.neo4j_label


with GraphDriver() as db:
    res = db.run(
        f"""
        MATCH (c1:{cmpd_label})
        RETURN c1 as chem
        """,
    ).value()


z={tuple([i[x] for x in ChemicalCompoundNode.required_properties]):ChemicalCompound(i['inchi']) for i in res}


big_vec=[i.morgan for i in z.values()]
cutoff=0.5

# loop through big_vec and compare each vector to the rest of the vectors
# if the similarity is greater than the cutoff, add the node to the list

big_dict = {i:[] for i in z.keys()}

id_list=list(z.keys())

chemsim_set = set()
for chemkey in id_list:
    for i, x in enumerate(
        DataStructs.BulkTanimotoSimilarity(z[chemkey].morgan,big_vec)):
        if x >= cutoff:
            if chemkey != id_list[i]:
                temp = [chemkey , id_list[i]]
                temp.sort()
                score = int(x * 100)
                chemsim_set.add(ChemicalSimilarity(start=z[temp[0]].node,end=z[temp[1]].node,properties={"similarity":score}))

ChemicalSimilarity.add_multiple_to_neo4j(list(chemsim_set),create=True)










bro = [
        (k for k in z.keys())
        for i, x in enumerate(
            DataStructs.BulkTanimotoSimilarity(big_vec[i],big_vec)
        )
        if x >= cutoff
    ]
