# # /home/chase/Downloads/rhea-chebi-smiles.tsv
# from collections import defaultdict
# from typing import List


# from socialgene.base.chem import ChemicalCompound
# from socialgene.neo4j.neo4j import GraphDriver
# from socialgene.utils.logging import log


# ChemicalFragments()

# from socialgene.base.chem import ChemicalCompound, ChemicalFragments

# a = []
# with open("/home/chase/Downloads/rhea-chebi-smiles.tsv") as f:
#     for l in f:
#         smiles = l.split("\t")[1]
#         a.append(ChemicalCompound(smiles))


# list(a.keys())[0]

# a["CHEBI:7"]
# b["CHEBI:7"]

# b = {k: Chem.MolToInchi(v.mol) for k, v in a.items()}


# bb = defaultdict(list)
# for k, v in b.items():
#     bb[v].append(k)

# bb = {k: v for k, v in bb.items() if len(v) > 1 and k != ""}


# result = set(
#     chain.from_iterable(values for key, values in rev_dict.items() if len(values) > 1)
# )


# with GraphDriver() as db:
#     results = db.run(
#         """
#         WITH $fragments as fragments
#         UNWIND fragments as fragment
#         MERGE (f:chemical_fragment {uid: fragment[0]})
#         MERGE (c:chemical_compound {uid: chemical_compound_uid})
#         MERGE (c)-[:CONTAINS]->(f)
#         """,
#         fragments=[(i, j) for i, j in a.fragments.items()],
#         chemical_compound_uid=a.uid,
#     ).value()


# with GraphDriver() as db:
#     results = db.run(
#         """
#         CREATE CONSTRAINT chemical_fragment
#         FOR (cf:chemical_compound) REQUIRE (cf.inchi, cf.CanonicalSmiles) IS NODE KEY
#         """,
#     ).value()


# def add_chemicals(props):
#     with GraphDriver() as db:
#         results = db.run(
#             """
#             WITH $props as props
#             UNWIND props as prop
#             MERGE (c:chemical_compound {inchi: prop.inchi, CanonicalSmiles: prop.CanonicalSmiles})
#             ON CREATE SET c = prop
#             """,
#             props=props,
#         ).value()


# add_chemicals([a.hash_dict | a.base_properties for a in a])


# with GraphDriver() as db:
#     results = db.run(
#         """
#         WITH $rows as rows
#         UNWIND rows as row
#         MATCH (c:chemical_compound {inchi: row[0], CanonicalSmiles: row[1]})
#         MERGE (f:chemical_fragment {uid: row[2]})
#         MERGE (c)-[:CONTAINS {n:row[3]}]->(f)
#         """,
#         rows=frag_list,
#     ).value()

# frag_list = []
# for i in a:
#     frag_list.extend(
#         [
#             (Chem.MolToInchi(i.mol), Chem.MolToSmiles(i.mol), k, v)
#             for k, v in i.fragments.items()
#         ]
#     )
