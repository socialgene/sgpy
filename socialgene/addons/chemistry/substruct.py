# from socialgene.neo4j.neo4j import GraphDriver
# import logging
# import os
# from socialgene.base.chem import ChemicalCompound
# from rich.progress import Progress
# from socialgene.addons.chemistry.nr import ChemicalSubstructure, ContainsSubstructure, ChemicalCompoundNode
# from rdkit import Chem
# from multiprocessing import Pool
# from rdkit.Chem.MolStandardize import rdMolStandardize
# from itertools import batched
# logging.getLogger("neo4j").setLevel(logging.WARNING)
# logging.getLogger().setLevel(logging.INFO)

# if __name__ == "__main__":

#     #inspect(a)
#     with GraphDriver() as db:
#         results = db.run(
#             """
#             MATCH (n:chemical_compound) RETURN n.inchi as inchi
#             """
#         ).value()


#     nodes=set()
#     rels=set()

#     def process_subgraph(batched_results):
#         nodes=set()
#         rels=set()
#         for i in batched_results:
#             sgmol = ChemicalCompound(i)
#             subgraphs = Chem.FindAllSubgraphsOfLengthN(sgmol.mol, 5)
#             for subgraph in subgraphs:
#                 sub_mol = Chem.PathToSubmol(sgmol.mol, subgraph, useQuery=True)
#                 sub_mol = rdMolStandardize.Cleanup(sub_mol)
#                 sub_mol = rdMolStandardize.Normalize(sub_mol)
#                 r = rdMolStandardize.Reionizer()
#                 sub_mol = r.reionize(sub_mol)
#                 Chem.RemoveStereochemistry( sub_mol )
#                 #sub_mol=Chem.MolToSmiles(sub_mol, canonical=True)
#                 node = ChemicalSubstructure()
#                 temp_compound = ChemicalCompound(sub_mol)
#                 node.fill_from_dict(temp_compound.base_properties | temp_compound.hash_dict)
#                 del temp_compound
#                 nodes.add(node)
#                 cn = ChemicalCompoundNode()
#                 cn.fill_from_dict(sgmol.base_properties | sgmol.hash_dict)
#                 rel = ContainsSubstructure(cn, node)
#                 rels.add(rel)
#         nodes[0].add_multiple_to_neo4j(nodes)
#         rels[0].add_multiple_to_neo4j(rels)

#     def update_progress(*args):
#         pg.update(task, advance=1)

#     with Progress(transient=True) as pg:
#         task = pg.add_task("Progress...", total=len(results))
#         batched_results = batched(results, 100)
#         with Pool() as p:
#             for i in p.imap_unordered(process_subgraph, batched_results):
#                 nodes.update(i[0])
#                 rels.update(i[1])
#                 update_progress()








