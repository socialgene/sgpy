# from socialgene.base.socialgene import SocialGene
# from socialgene.neo4j.neo4j import GraphDriver  # grab the the neo4j connection


# a = SocialGene()

# with GraphDriver() as db:
#     results = db.run(
#         """
# MATCH (p1:protein) RETURN p1.uid
# """
#     )
#     for i in results:
#         i = i[0]
#         _ = a.add_protein(uid=i)


# with GraphDriver() as db:
#     results = db.run(
#         """
# MATCH (a1:assembly) RETURN a1
# """
#     )
#     for i in results:
#         i = i[0]
#         a.add_assembly(id=i._properties["id"])
#         a.assemblies[i._properties["id"]].info = (
#             a.assemblies[i._properties["id"]].info | i._properties
#         )


# with GraphDriver() as db:
#     results = db.run(
#         """
# MATCH (n1:nucleotide)-[:ASSEMBLES_TO]-(a1:assembly) RETURN {n:n1, a:a1}
# """
#     )
#     for i in results:
#         i = i[0]
#         a.assemblies[i["a"]._properties["id"]].add_locus(
#             i["n"]._properties["internal_id"]
#         )
#         a.assemblies[i["a"]._properties["id"]].loci[
#             i["n"]._properties["internal_id"]
#         ].info = (
#             a.assemblies[i["a"]._properties["id"]]
#             .loci[i["n"]._properties["internal_id"]]
#             .info
#             | i["n"]._properties
#         )


# with GraphDriver() as db:
#     results = db.run(
#         """
# MATCH (a1:assembly)<-[:ASSEMBLES_TO]-(n1:nucleotide)-[e1:ENCODES]->(p1:protein) RETURN {a:a1.uid, n:n1.internal_id, p:p1.uid, start:e1.start, end:e1.end, strand:e1.strand}
# """
#     )
#     for i in results:
#         i = i[0]
#         a.assemblies[i["a"]].loci[i["n"]].add_feature(
#             type="protein",
#             id=i["p"],
#             start=i["start"],
#             end=i["end"],
#             strand=i["strand"],
#         )


# with GraphDriver() as db:
#     results = db.run(
#         """
# MATCH (p1:protein)<-[a1:ANNOTATES]-(h1:hmm) RETURN {p:p1.uid, a:a1, h:h1.uid}
# """
#     )
#     for i in results:
#         i = i[0]
#         z = {"hmm_id": i["h"]} | i["a"]._properties
#         a.proteins[i["p"]].add_domain(**z, exponentialized=False)

# with GraphDriver() as db:
#     results = db.run(
#         """
# MATCH (p1)-[b:BLASTP]->(p2:protein)
# RETURN
#     p1.uid as query,
#     p2.uid as subject,
#     b.bitscore as bitscore,
#     b.evalue as evalue,
#     b.gapopen as gapopen,
#     b.length as length,
#     b.mismatch as mismatch,
#     b.pident as pident,
#     b.qcovhsp as qcovhsp,
#     b.qend as qend,
#     b.qstart as qstart,
#     b.send as send,
#     b.sstart as sstart
# """
#     )
#     a.blastp_df = results.to_df()


# with GraphDriver() as db:
#     results = db.run(
#         """
# MATCH (p1)-[b:MMSEQS2]->(p2:protein)
# RETURN
#     p1.uid as cluster_representative,
#     p2.uid as external_id
# """
#     )
#     a.mmseqs2_df = results.to_df()


# with GraphDriver() as db:
#     results = db.run(
#         """
# MATCH (p1)-[b:SIMILAR]-(p2:protein)
# RETURN DISTINCT
#     p1.uid as p1,
#     p2.uid as p2,
#     b.score as score
# """
#     )
#     a.sim_df = results.to_df()
