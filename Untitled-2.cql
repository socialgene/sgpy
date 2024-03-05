






MATCH (bgc:nucleotide)-[:ENCODES]->(p1:protein)
WHERE bgc.external_id STARTS WITH "BGC"
WITH bgc, count(DISTINCT p1) as mm
SET bgc.distinct_proteins = mm



MATCH (bgc:nucleotide)
WHERE bgc.external_id STARTS WITH "BGC"
MATCH (bgc)-[:ENCODES]->(p0:protein)<-[:MMSEQS_90]-(:protein)<-[:MMSEQS_70]-(p1:protein)
MATCH (bgc2:nucleotide)-[:ENCODES]->(:protein)<-[:MMSEQS_90]-(:protein)<-[:MMSEQS_70]-(p1)
WITH bgc2, bgc, count(DISTINCT p0) as mm
WHERE mm / bgc.distinct_proteins   >= 0.7
CREATE (bgc)-[:SHARES_MMSEQS_70]->(bgc2)


MATCH (bgc:nucleotide)
WHERE bgc.external_id STARTS WITH "BGC"
MATCH (bgc)-[:ENCODES]->(p0:protein)<-[:MMSEQS_90]-(p1:protein)
MATCH (bgc2:nucleotide)-[e1:ENCODES]->(:protein)<-[:MMSEQS_90]-(p1)
WHERE bgc <> bgc2
WITH  bgc, bgc2, e1.antismash_region as region, count(DISTINCT p0) as matched_proteins
WHERE matched_proteins / bgc.distinct_proteins >= 0.7
MATCH (bgc)-[:ASSEMBLES_TO]-(a1:assembly)
MATCH (bgc2)-[:ASSEMBLES_TO]-(a2:assembly)
CREATE (bgc)-[:MIBIG_RELATIVES]->(A2)






WITH '
MATCH p=()-[r:MIBIG_RELATIVES]->() RETURN p
' AS query
CALL apoc.export.graphml.query(query, "import/movies-directed.graphml", {})
YIELD file, source, format, nodes, relationships, properties, time, rows, batchSize, batches, done, data
RETURN file, source, format, nodes, relationships, properties, time, rows, batchSize, batches, done, data;

