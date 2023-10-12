# import json
# from typing import List

# from socialgene.utils.logging import log


# def _flatten_list(x):
#     return sorted([item for row in x for item in row], reverse=True)


# class Clustermap:
#     def __init__(self, search_object):
#         self.uid = 0
#         self.uid_dict = {}
#         self.search_object = search_object

#     def _get_uid(self, obj=None):
#         # increment and return uid as string
#         self.uid += 1
#         self.uid_dict[str(self.uid)] = obj
#         return str(self.uid)

#     def _feature(self, feature_obj):
#         return {
#             "uid": self._get_uid(obj=feature_obj),
#             "label": feature_obj.protein_id,
#             "names": {
#                 "name": feature_obj.protein_id,
#                 "description": feature_obj.description,
#                 "hash": feature_obj.protein_hash,
#             },
#             "start": feature_obj.start,
#             "end": feature_obj.end,
#             "strand": feature_obj.strand,
#         }

#     def _locus(self, locus_name, locus_obj):
#         return {
#             "uid": self._get_uid(obj=locus_obj),
#             "name": locus_name,
#             "genes": [self._feature(feature_obj) for feature_obj in locus_obj.features],
#             "start": min([i.start for i in locus_obj.features]),
#             "end": max([i.end for i in locus_obj.features]),
#         }

#     def _loci(self, assembly_obj):
#         return [
#             self._locus(locus_name=k, locus_obj=v) for k, v in assembly_obj.loci.items()
#         ]

#     def _clusters(self, bgc_order):
#         log.info("Creating clustermap.js clusters")
#         return {
#             "clusters": [
#                 {
#                     "uid": self._get_uid(obj=self.sg_object.assemblies[k]),
#                     "name": self.sg_object.assemblies[k].name,
#                     "loci": self._loci(self.sg_object.assemblies[k]),
#                 }
#                 for k in bgc_order
#             ]
#         }

#     def _get_gene_uid_of_protein_hash(self, protein_hash: str) -> List[str]:
#         """Search the self.uid_dict for values that have a given protein hash

#         Args:
#             protein_hash (str): protein uid

#         Returns:
#             List[str]: _description_
#         """
#         return [
#             k
#             for k, v in self.uid_dict.items()
#             if type(v).__name__ == "Feature" and v.protein_hash == protein_hash
#         ]

#     def _create_group(self, k, v):
#         """Create a clustermap group"""
#         return {
#             "uid": self._get_uid(obj=k),
#             "label": group_dict_info[k][1],
#             "genes": _flatten_list([self._get_gene_uid_of_protein_hash(i) for i in v]),
#         }

#     def _create_groups_json(self, groupdict, group_dict_info):
#         log.info("Creating clustermap.js groups")
#         return {
#             "groups": [
#                 self._create_group(group_dict_info, k, v) for k, v in groupdict.items()
#             ]
#         }

#     def _create_links(self, df):
#         log.info("Creating clustermap.js links")
#         res = list()
#         # ugly because it first filters the protein_comparison table for proteins present in the self.sg_object
#         for index, row in df.iterrows():
#             for query_uid in self._get_gene_uid_of_protein_hash(row["query"]):
#                 for target_uid in self._get_gene_uid_of_protein_hash(row.target):
#                     # don't do self-links they mess up clustermap.js groups
#                     if (
#                         self.uid_dict[query_uid].parent_object
#                         != self.uid_dict[target_uid].parent_object
#                     ):
#                         res.append(
#                             {
#                                 "uid": self._get_uid(),
#                                 "query": {
#                                     "uid": query_uid,
#                                     "name": self._get_uid(),
#                                 },
#                                 "target": {
#                                     "uid": target_uid,
#                                     "name": self._get_uid(),
#                                 },
#                                 "identity": round(row.score * 2 / 3, 2),
#                             }
#                             # row["mod_score"] * 2/3 adjusts max mod score of 1.5 to ~1
#                         )
#         return res

#     def _build(self, df, bgc_order, group_df):
#         return (
#             self._clusters(bgc_order)
#             | self._create_groups_json(group_df)
#             | {"links": self._create_links(df)}
#         )

#     def write(self, bgc_df, bgc_order, link_df, group_df, outpath):
#         log.info(f"Writing clustermap.js output to: {outpath}")
#         with open(outpath, "w") as outfile:
#             json.dump(
#                 self._build(
#                     bgc_order=bgc_order,
#                     link_df=link_df,
#                     group_df=group_df,
#                 ),
#                 outfile,
#             )
