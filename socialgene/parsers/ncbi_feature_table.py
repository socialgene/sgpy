#
# import gzip
# from mimetypes import guess_type
# from functools import partial
# import sqlite3

#

#
# from socialgene.utils.logging import log
# from socialgene.base.socialgene import SocialGene
# from sqlite3 import Error

# # see assumption note in ./socialgene/src/socialgene/utils/protein_sqlite.py


# class NCBIFeatureTable(SocialGene):
#     def __init__(self):
#         super().__init__()
#         self.conn = None

#     def connect(self, db_path):
#         try:
#             self.conn = sqlite3.connect(db_path)
#         except Error as e:
#             log.error(e)

#     def sql_hash_lookup(self, accession):
#         cur = self.conn.cursor()
#         cur.execute("SELECT sha512t24u FROM protmap WHERE accession=?", ([accession]))
#         rows = cur.fetchall()
#         # ensure only one result was returned
#         if not rows:
#             raise ValueError(
#                 f'Couldn\'t find a match for "{accession}" in the sqlite database'
#             )
#         if not len(rows[0]) == 1:
#             raise ValueError(
#                 f'Found multiple hashes for "{accession}" in the sqlite database'
#             )
#         return rows[0][0]

#     def parse(self, input_path):
#         encoding = guess_type(input_path)[1]
#         if encoding == "gzip":
#             _open = partial(gzip.open, mode="rt")
#         else:
#             _open = open
#         with _open(input_path) as f:
#             for line in f:
#                 if line.startswith("#"):
#                     continue
#                 temp = line.split("\t")
#                 self.add_assembly(id=temp[2])
#                 self.assemblies[temp[2]].add_locus(id=temp[6])
#                 uid = self.sql_hash_lookup(temp[10])
#                 self.assemblies[temp[2]].loci[temp[6]].add_feature(
#                     type="protein",
#                     id=uid,
#                     start=temp[7],
#                     end=temp[8],
#                     strand=temp[9],
#                 )
#                 self.add_protein(
#                     description=temp[13],
#                     external_id=temp[10],
#                     uid=uid,
#                 )
#         log.info(f"Parsed {len(self.assemblies)} assemblies")
#         log.info(f"Parsed {sum([len(i.loci) for i in self.assemblies.values()])} loci")
#         log.info(f"Parsed {len(self.proteins)} non-redundant proteins")


# # accession="BGC0000001|c1|1-1083|+|AEK75490.1|protein_methyltransferase|AEK75490.1"
