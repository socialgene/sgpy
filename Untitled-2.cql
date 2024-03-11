






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
MATCH (bgc)-[:ENCODES]->(p0:protein)<-[:MMSEQS_90]-(:protein)<-[:MMSEQS_70]-(p1:protein)
MATCH (bgc2:nucleotide)-[e1:ENCODES]->(:protein)<-[:MMSEQS_90]-(:protein)<-[:MMSEQS_70]-(p1)
WHERE bgc <> bgc2
WITH  bgc, bgc2, e1.antismash_region as region, count(DISTINCT p0) as matched_proteins
WHERE matched_proteins / bgc.distinct_proteins >= 0.7
MATCH (bgc)-[:ASSEMBLES_TO]-(a1:assembly)
MATCH (bgc2)-[:ASSEMBLES_TO]-(a2:assembly)
WITH bgc, bgc2, toFloat(matched_proteins) / toFloat(bgc.distinct_proteins) as score
CREATE (a1)-[:MIBIG_RELATIVES {score:score}]->(a2)








:auto MATCH (bgc:nucleotide)
WHERE bgc.external_id STARTS WITH "BGC"
CALL {
  WITH bgc
  MATCH (bgc)-[:ENCODES]->(p0:protein)<-[:MMSEQS_90]-(:protein)<-[:MMSEQS_70]-(p1:protein)
  MATCH (bgc2:nucleotide)-[e1:ENCODES]->(:protein)<-[:MMSEQS_90]-(:protein)<-[:MMSEQS_70]-(p1)
  WHERE bgc <> bgc2 and not bgc2.external_id STARTS WITH "BGC"
  WITH  bgc, bgc2, e1.antismash_region as region, count(DISTINCT p0) as matched_proteins
  WITH bgc, bgc2, toFloat(matched_proteins) / toFloat(bgc.distinct_proteins) as score
  WHERE score >= 0.7
  MATCH (bgc)-[:ASSEMBLES_TO]->(a1:assembly)
  MATCH (bgc2)-[:ASSEMBLES_TO]->(a2:assembly)
  CREATE (a1)-[:MIBIG_ANTISMASH_RELATIVES {score:score}]->(a2)

} IN TRANSACTIONS OF 100 ROWS
RETURN 1;



:auto MATCH (bgc:nucleotide)
WHERE bgc.external_id STARTS WITH "BGC"
CALL {
  WITH bgc
  MATCH (bgc)-[:ENCODES]->(p0:protein)<-[:MMSEQS_90]-(:protein)<-[:MMSEQS_70]-(:protein)<-[:MMSEQS_50]-(p1:protein)
  MATCH (bgc2:nucleotide)-[e1:ENCODES]->(:protein)<-[:MMSEQS_90]-(:protein)<-[:MMSEQS_70]-(:protein)<-[:MMSEQS_50]-(p1)
  WHERE bgc <> bgc2 and not bgc2.external_id STARTS WITH "BGC"
  WITH  bgc, bgc2, e1.antismash_region as region, count(DISTINCT p0) as matched_proteins
  WITH bgc, bgc2, toFloat(matched_proteins) / toFloat(bgc.distinct_proteins) as score
  WHERE score >= 0.7
  MATCH (bgc)-[:ASSEMBLES_TO]->(a1:assembly)
  MATCH (bgc2)-[:ASSEMBLES_TO]->(a2:assembly)
  CREATE (a1)-[:MIBIG_ANTISMASH_RELATIVES_50 {score:score}]->(a2)

} IN TRANSACTIONS OF 50 ROWS
RETURN 1;



:auto MATCH (bgc:nucleotide)
WHERE bgc.external_id STARTS WITH "BGC"
CALL {
  WITH bgc
  MATCH (bgc)-[:ENCODES]->(p0:protein)<-[:MMSEQS_90]-(:protein)<-[:MMSEQS_70]-(p1:protein)
  MATCH (a2:assembly)<-[:ASSEMBLES_TO]-(:nucleotide)-[e1:ENCODES]->(:protein)<-[:MMSEQS_90]-(:protein)<-[:MMSEQS_70]-(p1)
  WHERE not a2.uid STARTS WITH "BGC"
  WITH  bgc, a2, count(DISTINCT p0) as matched_proteins
  WITH bgc, a2, toFloat(matched_proteins) / toFloat(bgc.distinct_proteins) as score
  WHERE score >= 0.7
  MATCH (bgc)-[:ASSEMBLES_TO]->(a1:assembly)
  CREATE (a1)-[:MIBIG_RELATIVES {score:score}]->(a2)
} IN TRANSACTIONS OF 100 ROWS
RETURN 1;





:auto MATCH (bgc:nucleotide)
WHERE bgc.external_id STARTS WITH "BGC"
CALL {
  WITH bgc
  MATCH (bgc)-[:ENCODES]->(p0:protein)<-[:MMSEQS_90]-(:protein)<-[:MMSEQS_70]-(:protein)<-[:MMSEQS_50]-(p1:protein)
  MATCH (a2:assembly)<-[:ASSEMBLES_TO]-(:nucleotide)-[e1:ENCODES]->(:protein)<-[:MMSEQS_90]-(:protein)<-[:MMSEQS_70]-(:protein)<-[:MMSEQS_50]-(p1)
  WHERE not a2.uid STARTS WITH "BGC"
  WITH  bgc, a2, count(DISTINCT p0) as matched_proteins
  WITH bgc, a2, toFloat(matched_proteins) / toFloat(bgc.distinct_proteins) as score
  WHERE score >= 0.7
  MATCH (bgc)-[:ASSEMBLES_TO]->(a1:assembly)
  CREATE (a1)-[:MIBIG_RELATIVES_50 {score:score}]->(a2)
} IN TRANSACTIONS OF 50 ROWS
RETURN 1;





WITH '
MATCH p=()-[r:MIBIG_RELATIVES]->() RETURN p
' AS query
CALL apoc.export.graphml.query(query, "import/MIBIG_RELATIVES.graphml", {})
YIELD file, source, format, nodes, relationships, properties, time, rows, batchSize, batches, done, data
RETURN file, source, format, nodes, relationships, properties, time, rows, batchSize, batches, done, data;

WITH '
MATCH p=()-[r:MIBIG_RELATIVES_50]->() RETURN p
' AS query
CALL apoc.export.graphml.query(query, "import/MIBIG_RELATIVES_50.graphml", {})
YIELD file, source, format, nodes, relationships, properties, time, rows, batchSize, batches, done, data
RETURN file, source, format, nodes, relationships, properties, time, rows, batchSize, batches, done, data;

WITH '
MATCH p=()-[r:MIBIG_ANTISMASH_RELATIVES]->() RETURN p
' AS query
CALL apoc.export.graphml.query(query, "import/MIBIG_ANTISMASH_RELATIVES.graphml", {})
YIELD file, source, format, nodes, relationships, properties, time, rows, batchSize, batches, done, data
RETURN file, source, format, nodes, relationships, properties, time, rows, batchSize, batches, done, data;

WITH '
MATCH p=()-[r:MIBIG_ANTISMASH_RELATIVES_50]->() RETURN p
' AS query
CALL apoc.export.graphml.query(query, "import/MIBIG_ANTISMASH_RELATIVES_50.graphml", {})
YIELD file, source, format, nodes, relationships, properties, time, rows, batchSize, batches, done, data
RETURN file, source, format, nodes, relationships, properties, time, rows, batchSize, batches, done, data;




Set 383812 properties, created 383812 relationships, started streaming 2502 records after 36 ms and completed after 163173 ms, displaying first 1000 rows.
