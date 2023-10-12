# import json
# from socialgene.base.socialgene import SocialGene
# from socialgene.clustermap.clustermap import Clustermap
# import os
# import tempfile
# from pathlib import Path


# sg_object = SocialGene()
# sg_object.add_assembly("a")
# ######################################################
# sg_object.assemblies["a"].add_locus("a")
# sg_object.assemblies["a"].loci["a"].add_feature(
#     protein_hash="p1-a",
#     type="protein",
#     start=1,
#     end=1000,
#     strand=1,
# )
# sg_object.assemblies["a"].loci["a"].add_feature(
#     protein_hash="p2",
#     type="protein",
#     start=1100,
#     end=2000,
#     strand=-1,
# )
# ######################################################
# sg_object.add_assembly("b")
# sg_object.assemblies["b"].add_locus("b")
# sg_object.assemblies["b"].loci["b"].add_feature(
#     protein_hash="p1-b",
#     type="protein",
#     start=1,
#     end=1000,
#     strand=1,
# )
# sg_object.assemblies["b"].loci["b"].add_feature(
#     protein_hash="p2",
#     type="protein",
#     start=1100,
#     end=2000,
#     strand=-1,
# )
# ######################################################
# sg_object.add_assembly("c")
# sg_object.assemblies["c"].add_locus("c")
# sg_object.assemblies["c"].loci["c"].add_feature(
#     protein_hash="p1-c",
#     type="protein",
#     start=1,
#     end=1000,
#     strand=1,
# )
# sg_object.assemblies["c"].loci["c"].add_feature(
#     protein_hash="p2",
#     type="protein",
#     start=1100,
#     end=2000,
#     strand=-1,
# )


# cmap = Clustermap()
# groupdict = {
#     "g1": [
#         "p1-a",
#         "p1-b",
#         "p1-c",
#     ],
#     "g2": [
#         "p2",
#     ],
# }
# group_dict_info = {
#     "g1": ("g1", "group1"),
#     "g2": ("g1", "group2"),
# }

# assembly_order = [
#     "a",
#     "c",
#     "b",
# ]

# ## Which proteins match to the input, what happens if match to more than 1?


# #   IDS ARE BEING ASSIGNED DIFFERENTLY EACH RUN


# def test_clustermap_write():
#     with tempfile.NamedTemporaryFile() as fp:
#         cmap.write(
#             sg_object,
#             groupdict=groupdict,
#             group_dict_info=group_dict_info,
#             assembly_order=assembly_order,
#             outpath=fp,
#         )
#         with open(fp.name, "r") as h:
#             data = json.load(h)
#     assert data == {
#         "clusters": [
#             {
#                 "uid": "1",
#                 "name": "a",
#                 "loci": [
#                     {
#                         "uid": "2",
#                         "name": "a",
#                         "genes": [
#                             {
#                                 "uid": "3",
#                                 "label": None,
#                                 "names": {
#                                     "name": None,
#                                     "description": None,
#                                     "hash": "p2",
#                                 },
#                                 "start": 1100,
#                                 "end": 2000,
#                                 "strand": -1,
#                             },
#                             {
#                                 "uid": "4",
#                                 "label": None,
#                                 "names": {
#                                     "name": None,
#                                     "description": None,
#                                     "hash": "p1-a",
#                                 },
#                                 "start": 1,
#                                 "end": 1000,
#                                 "strand": 1,
#                             },
#                         ],
#                         "start": 1,
#                         "end": 2000,
#                     }
#                 ],
#             },
#             {
#                 "uid": "5",
#                 "name": "c",
#                 "loci": [
#                     {
#                         "uid": "6",
#                         "name": "c",
#                         "genes": [
#                             {
#                                 "uid": "7",
#                                 "label": None,
#                                 "names": {
#                                     "name": None,
#                                     "description": None,
#                                     "hash": "p1-c",
#                                 },
#                                 "start": 1,
#                                 "end": 1000,
#                                 "strand": 1,
#                             },
#                             {
#                                 "uid": "8",
#                                 "label": None,
#                                 "names": {
#                                     "name": None,
#                                     "description": None,
#                                     "hash": "p2",
#                                 },
#                                 "start": 1100,
#                                 "end": 2000,
#                                 "strand": -1,
#                             },
#                         ],
#                         "start": 1,
#                         "end": 2000,
#                     }
#                 ],
#             },
#             {
#                 "uid": "9",
#                 "name": "b",
#                 "loci": [
#                     {
#                         "uid": "10",
#                         "name": "b",
#                         "genes": [
#                             {
#                                 "uid": "11",
#                                 "label": None,
#                                 "names": {
#                                     "name": None,
#                                     "description": None,
#                                     "hash": "p1-b",
#                                 },
#                                 "start": 1,
#                                 "end": 1000,
#                                 "strand": 1,
#                             },
#                             {
#                                 "uid": "12",
#                                 "label": None,
#                                 "names": {
#                                     "name": None,
#                                     "description": None,
#                                     "hash": "p2",
#                                 },
#                                 "start": 1100,
#                                 "end": 2000,
#                                 "strand": -1,
#                             },
#                         ],
#                         "start": 1,
#                         "end": 2000,
#                     }
#                 ],
#             },
#         ],
#         "groups": [
#             {"uid": "13", "label": "group1", "genes": ["4", "11", "7"]},
#             {"uid": "14", "label": "group2", "genes": ["3", "8", "12"]},
#         ],
#         "links": [],
#     }
