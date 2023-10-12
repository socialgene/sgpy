// Name: minimal_example
// Description:
// Param: ["WP_139587022.1", "WP_139587010.1"]
WITH $param AS input_param
MATCH (n:protein)
WHERE n.name IN input_param
RETURN n.uid
LIMIT 25;

// Name: total_node_count
// Description: Count all nodes
// Param:
MATCH (n)
RETURN count(n) AS count;

// Name: total_relationship_count
// Description: Count all relationships
// Param:
MATCH ()-[r]->()
RETURN count(r) AS count;

// Name: node_label_count
// Description: Get count for each node label
// Param:
CALL db.labels() YIELD label
CALL apoc.cypher.run('MATCH (:`'+label+'`) RETURN count(*) as count', { }) YIELD value
RETURN apoc.map.fromPairs(collect([label, value.count])) AS count;

// Name: relationship_label_count
// Description: Get count for each relationship label
// Param:
CALL db.relationshipTypes() YIELD relationshipType AS type
CALL apoc.cypher.run('MATCH ()-[:`'+type+'`]->() RETURN count(*) as count', { }) YIELD value
RETURN apoc.map.fromPairs(collect([type, value.count])) AS count;

// Name: database_parameters
// Description: Get info about the parameters used when creating the database
// Param:
MATCH (p:parameters)
RETURN p;

// Name: search_protein_hash
// Description:
// Param: ["8Q02h-w20QNX9aFt1mx2mZgGcXyQGCMn"]
OPTIONAL MATCH (n:protein)
WHERE n.uid IN $param
RETURN collect(n.uid) AS result;

// Name: get_assembly_ids
// Description:
// Param:
MATCH (a1:assembly)
RETURN collect(a1.uid) AS assemblies;

// Name: find_identical_domain_content
// Description: Finds a single protein with identical domain content for an input list of domains
// Param:
WITH $param AS input_domains
MATCH (h1:hmm)-[:ANNOTATES]->(p1:protein)<-[:ANNOTATES]-(h2:hmm)
WHERE h1.uid IN input_domains
WITH p1, collect( DISTINCT h2.uid) AS domains_per_protein, input_domains
WITH p1, apoc.coll.isEqualCollection(domains_per_protein, input_domains) AS is_equal, domains_per_protein, input_domains
RETURN p1.uid AS protein_hash
 ORDER BY is_equal DESC
LIMIT 1;

// Name: find_similar_bgc
// Description: Finds similar bgc as input protein domains
// Param: [{prot:'', ids:[]}]
WITH $param AS input_proteins
WITH input_proteins, size(input_proteins) AS size_input_proteins
UNWIND input_proteins AS single_protein
call{
  WITH single_protein, size_input_proteins
  MATCH (h0:hmm)
  WHERE h0.uid IN single_protein['domains']
  WITH single_protein, single_protein['domains'][0..10] AS input_subset, COLLECT(id(h0)) AS input_domains_db_id
  MATCH (prot1:protein)<-[a1:ANNOTATES]-(hmm1:hmm)
  WHERE hmm1.uid IN input_subset
  WITH single_protein, COLLECT( DISTINCT (hmm1.uid)) AS subset_match, input_domains_db_id, prot1
  WHERE size(subset_match) > size(single_protein['domains']) *.8 AND size(subset_match) < size(single_protein['domains']) * 1.2
  MATCH (prot1)<-[a2:ANNOTATES]-(hmm2:hmm)
  WITH single_protein, prot1, COLLECT( DISTINCT (id(hmm2))) AS tmp, input_domains_db_id
  WHERE gds.similarity.jaccard(tmp, input_domains_db_id) > .80
  MATCH (hmms3:hmm)-[:ANNOTATES]-(prot1)<-[c1:ENCODES]-(n1:nucleotide)-[:ASSEMBLES_TO]->(a1:assembly)
  RETURN hmms3 AS result_hmms, a1.uid AS assembly, prot1 AS match_protein, n1.uid AS locus, c1.start AS locus_start, c1.end AS locus_end, c1.strand AS strand, single_protein['prot'] AS query_protein
}
WITH assembly, { query_protein:query_protein, locus:locus, match_protein_id:match_protein.uid, match_protein_name:match_protein.name, locus_end:locus_end, locus_start:locus_start, strand:strand } AS a1, query_protein, size_input_proteins
WITH assembly, collect( DISTINCT (query_protein)) AS distinct_query_protein, collect(a1) AS locus_info, size_input_proteins
RETURN assembly, distinct_query_protein, locus_info;

// Name: get_species_and_assembly_from_protein_list
// Description:
// Param:
WITH $param AS input_protein_list
MATCH (p1:protein)<-[:ENCODES]-(:nucleotide)-[:ASSEMBLES_TO]->(a1:assembly)-[:IS_TAXON]-(t1:taxid)
WHERE p1.uid IN input_protein_list
RETURN p1.uid AS protein, collect( DISTINCT (a1.uid)) AS assemblies, collect( DISTINCT (t1.name)) AS species;

// Name: get_class_taxid_from_assembly_id
// Description:
// Param:
WITH $param AS input_list
MATCH (n:assembly)-[:IS_TAXON]-(:taxid)-[:TAXON_PARENT*1..6]->(t1:taxid)
WHERE n.uid IN input_list AND t1.rank = "class"
RETURN n.uid AS assembly_id, t1.name AS class_taxon

// Name: get_genus_taxid_from_assembly_id
// Description:
// Param:
WITH $param AS input_list
MATCH (n:assembly)-[:IS_TAXON]-(:taxid)-[:TAXON_PARENT*1..4]->(t1:taxid)
WHERE n.uid IN input_list AND t1.rank = "genus"
RETURN n.uid AS assembly_id, t1.name AS class_taxon

// Name: get_protein_name
// Description: Get the original protein accessions, provided a list of protein hashes
// Param:
WITH $param AS input_protein_list
MATCH (p1:protein)
WHERE p1.uid IN input_protein_list
RETURN p1.uid AS protein_id, p1.name AS target_id, p1.description AS description;

// Name: get_protein_info
// Description: Get the original protein accessions, provided a list of protein hashes
// Param:
WITH $param AS input_protein_list
MATCH (p1:protein)
WHERE p1.uid IN input_protein_list
RETURN p1.uid AS protein_id, p1.name AS name, p1.description AS description;

// Name: get_blastp_matches
// Description:
// Param:
WITH $param AS input
MATCH (p1:protein)
WHERE p1.uid IN input
MATCH (p1)-[b1:BLASTP]->(p2:protein)
RETURN p1.uid, p2.uid, p1.name, p2.name, p1.description, p2.description, b1.`:bitscore` AS blast_bitscore, b1.`:pident` AS p_ident, b1.`:evalue` AS blast_evalue, b1.`:length` AS blast_length, p1.seqlen AS p1_seqlen, p2.seqlen AS p2_seqlen;

// Name: get_protein_domains
// Description:
// Param:
WITH $param AS input
MATCH (p1:protein)-[a1:ANNOTATES]-(h1:hmm)
WHERE p1.uid IN input
RETURN p1.uid, collect({ hmm_id:h1.uid, domain_properties:properties(a1) }) AS domains;

// Name: get_hmms_given_protein_ids
// Description:
// Param: list of hmms
MATCH (prot1:protein)<-[a1:ANNOTATES]-(hmm1:hmm)
WHERE prot1.uid IN $param
WITH prot1, collect(properties(a1)) AS hitlist, hmm1
RETURN prot1.uid AS protein, { hmms:collect({hmm_id:hmm1.uid, hitlist:hitlist })} AS domains;

// Name: get_protein_from_assembly
// Description:
// Param:
WITH $param AS input
MATCH (p1:protein)<-[:ENCODES]-(:nucleotide)-[:ASSEMBLES_TO]->(a1:assembly)
WHERE a1.uid IN input
RETURN collect(p1.uid) AS proteins;

// Name: get_assembly_from_proteins
// Description:
// Param:
WITH $param AS input
MATCH (p1:protein)<-[:ENCODES]-(:nucleotide)-[:ASSEMBLES_TO]->(a1:assembly)
WHERE p1.uid IN input
RETURN collect( DISTINCT (a1.uid)) AS assemblies;

// Name: get_assembly_and_locus_from_protein
// Description:
// Param:
WITH $param AS input
MATCH (p1:protein)<-[:ENCODES]-(n1:nucleotide)-[:ASSEMBLES_TO]->(a1:assembly)
WHERE p1.uid IN input
WITH a1, collect( DISTINCT (n1.uid)) AS loci
RETURN a1.uid AS assembly, loci;

// Name: count_matched_query_proteins_per_assembly
// Description:
// Param:
WITH $param AS inputs
UNWIND inputs AS input
MATCH (p1:protein)<-[:ENCODES]-(n1:nucleotide)-[:ASSEMBLES_TO]->(a1:assembly)
WHERE p1.uid IN input['target']
WITH a1, n1, apoc.coll.frequenciesAsMap(collect(input['query'])) AS counts
RETURN a1.uid, collect(n1.uid), counts

// Name: get_all_protein_ids
// Description:
// Param: None
MATCH (p1:protein)
RETURN collect(p1.uid) AS protein_hash_ids;

// Name: get_tuples_of_proteins_with_shared_domains
// Description:
// Param: None
MATCH (p1:protein)-[a1:ANNOTATES]-(h1:hmm)-[a2:ANNOTATES]-(p2:protein)
WHERE p1 <> p2
WITH collect([p1.uid, p2.uid]) AS tupe
RETURN tupe;

// Name: get_hmm_hashes_for_protein
// Description:
// Param:
WITH $param AS z
MATCH (p:protein)<-[r1:ANNOTATES]-(h:hmm)
WHERE p.uid IN z AND r1.ievalue < .1 AND r1.score1 > 20
WITH p, r1, h
 ORDER BY (r1.from_seq + r1.to_seq) / 2 //order by the midpoint of annotation
WITH p, collect(h.uid) AS hmms
RETURN (p.uid) AS protein_hash, hmms AS hmm_hashes;

// Name: get_bgcs
// Description:
// Param:
WITH $param AS a
MATCH (p1:protein)-[c1:ENCODES]-(n1:nucleotide)-[:ASSEMBLES_TO]-(a1:assembly)
WHERE p1.uid IN a
WITH a1, n1, c1, { id: p1.name, start:c1.start, end:c1.end, strand:c1.strand } AS tmp2
WITH a1, { uid:n1.uid, name:n1.uid, genes:collect(tmp2), start:min(c1.start), end:max(c1.end) } AS tmp4
RETURN { uid:ID(a1), name:a1.uid, loci:collect(tmp4) } AS clusters;

// Name: get_hmms_for_domain_class
// Description:
// Param:
WITH $param AS z
MATCH (p:protein)<-[r1:ANNOTATES]-(h:hmm)
WHERE p.uid IN z //AND r1.ievalue < .1 AND r1.score1 > 20
WITH p, r1, h
 ORDER BY (r1.from_seq + r1.to_seq) / 2 //order by the midpoint of annotation
WITH p, collect({ hash_id:h.uid, ali_from:r1.from_seq, ali_to:r1.to_seq, i_evalue:r1.ievalue }) AS domain_info
RETURN (p.uid) AS protein_hash, domain_info AS domain_info;

// Name: find_clusterjs_as_pandas_identicals
// Description:
// Param:
WITH $param AS input_groups
MATCH (p1:protein)-[c1:ENCODES]-(n1:nucleotide)-[:ASSEMBLES_TO]-(a1:assembly)
WHERE p1.uid IN input_groups
RETURN a1.uid AS assembly, n1.uid AS locus, p1.uid AS gene_hash, p1.name AS gene_name, c1.start AS start, c1.end AS end, (c1.start + c1.end) / 2 AS median_position, p1.uid AS
MATCH ;

// Name: query_1
// Description:
// Param:
WITH $param AS z
WITH z, size(z) AS y
UNWIND z AS x
call {
  WITH x
  MATCH (h1:hmm)-[r1:ANNOTATES]->(p1:protein)<-[c1:ENCODES]-(n1:nucleotide)-[:ASSEMBLES_TO]->(a1:assembly)
  WHERE id(h1) IN x['ids'] AND r1.ievalue < .1 AND r1.score1 > 2
  WITH collect(id(h1)) AS tmp, p1, a1, x, c1.end AS endc, c1.start AS startc, c1.strand AS strand , n1.uid AS nuc, x['prot'] AS query_protein
  WHERE size(tmp) > (0.8 * x['length']) AND size(tmp) < (1.2 * x['length'])
  RETURN a1.uid AS assembly, p1, nuc, startc, endc, strand, query_protein
}
call {
  WITH x, assembly, p1, nuc, startc, endc, strand, query_protein
  MATCH (h2:hmm)-[r2:ANNOTATES]->(p1)
  WHERE r2.ievalue < .1 AND r2.score1 > 2
  WITH collect(id(h2)) AS result_hmms, p1, assembly, nuc, startc, endc, strand, query_protein, x
  WHERE size(result_hmms) > (0.8 * x['length']) AND size(result_hmms) < (1.2 * x['length'])
  RETURN p1.name AS match_protein, endc AS locus_end, startc AS locus_start, nuc AS locus, result_hmms AS result_hmms
}
WITH assembly, collect([query_protein, locus, match_protein, locus_end, locus_start, strand, result_hmms]) AS a1, collect( DISTINCT query_protein) AS a2, y
WHERE size(a2) > 0.5 * y
RETURN assembly, a1;

// Name: find_similar_bgc3
// Description:
// Param: [{prot:'', ids:[]}]
WITH $param AS input_proteins
UNWIND input_proteins AS single_protein
WITH single_protein['prot'] AS input_protein_hash, single_protein['domains'] AS input_protein_domains
MATCH (h0:hmm)
WHERE h0.uid IN input_protein_domains
MATCH (prot1:protein)<-[a1:ANNOTATES]-(h0)
WITH input_protein_hash, size(input_protein_domains) AS tmp, prot1, size(collect( DISTINCT (id(h0)))) AS hmm_matches
WHERE abs(tmp-hmm_matches) < tmp * 0.1
MATCH (prot1)<-[a2:ANNOTATES]-(hmm2:hmm)
WITH input_protein_hash, tmp, prot1, size(collect( DISTINCT (id(hmm2)))) AS hmm_matches2
WHERE abs(tmp-hmm_matches2) < tmp * 0.1
RETURN input_protein_hash AS query, collect( DISTINCT (prot1.uid)) AS target;

// Name: search_a_single_protein
// Description:
// Param: [{prot:'', ids:[]}]
WITH $param AS input_protein_domains
WITH input_protein_domains, size(input_protein_domains) AS input_len
MATCH (h1:hmm)-[a2:ANNOTATES]->(prot1:protein)<-[a1:ANNOTATES]-(h0:hmm)
WHERE h0.uid IN input_protein_domains
WITH input_len, prot1, count(DISTINCT (h0)) AS hmm_matches,  count(DISTINCT (h1)) AS all_annotations
WHERE abs(input_len-hmm_matches) < 1 AND abs(input_len-all_annotations) < 1
RETURN collect(DISTINCT prot1.uid) AS target

// Name: retrieve_protein_locations
// Description: Given a list of proteins, retrieve all loci and assemblies that contain them
// Param: []
WITH $param AS input_protein_ids
MATCH (prot1:protein)<-[c1:ENCODES]-(n1:nucleotide)-[:ASSEMBLES_TO]->(a1:assembly)
WHERE prot1.uid IN input_protein_ids
WITH a1, n1, collect({ protein_id:prot1.uid, locus_start:c1.start, locus_end:c1.end, strand:c1.strand }) AS tmp
RETURN a1.uid AS assembly, collect({ locus:n1.uid, features:tmp }) AS loci;

// Name: temp1
// Description:
// Param: []
MATCH (p1)-[r:SIMILAR_jaccard]->(p2)
WITH r, [p1.uid, p2.uid] AS tmp
RETURN collect(tmp) AS res;

// Name: temp2
// Description:
// Param: []
WITH $param AS input
MATCH (p1), (p2)
WHERE p1.uid = input[0] AND p2.uid = input[1]
MERGE (p1)-[:SCORES { l1: input[2], l2: input[3], lev: input[4], jacc: input[5], mod: input[6] }]-(p2)

// Name: get_mmseqs_cluster_members
// Description: Retreive all proteins connected to an MMseqs2 cluster representative, self hits are removed
// Param: []
WITH $param AS input
MATCH (p1:protein)-[:MMSEQS2]->(p2:protein)
WHERE p1.uid IN input AND p1 <> p2
RETURN p1.uid AS cluster_representative, p2.uid AS member;





