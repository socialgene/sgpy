# a

1. Take an input BGC
2. Parse all the proteins in the input BGC
3. Annotate every input BGC proteins with HMM models, either pulling from Neo4j if the exact protein exists, or with HMMER
4. Get the list of all HMM models annotating input BGC proteins
5. Query the database to find the outdegree (calculate/populate if necessary) for each HMM model. (:hmm)-[:ANNOTATES]->(:protein)
6. Prioritize whcih input proteins to search based on the outdegree of its HMM annotations
    - max_query_proteins (float): Max proteins to return. If >0 and <1, will return X% of input proteins
    - max_domains_per_protein (int): Max domains to retain for each individual protein (highest outdegree dropped first)
    - max_outdegree (int): HMM model annotations with an outdegree higher than this will be dropped
    - scatter (bool, optional): Choose a random subset of proteins to search that are spread across the length of the input BGC. Defaults to False.
    - bypass (List[str], optional): List of locus tags that will bypass filtering. This is the ID found in a GenBank file "/locus_tag=" field. Defaults to None.
    - protein_id_bypass_list (List[str], optional): Less preferred than `bypass`. List of external protein IDs that will bypass filtering. This is the ID found in a GenBank file "/protein_id=" field. Defaults to None.
7. Search the database for all proteins that have the same HMM model annotations as the input BGC proteins
    - Output from database is a data frame with columns: ['assembly_uid', 'nucleotide_uid', 'target', 'n_start', 'n_end', 'query']
8. The initial hits output is filtered based on the following criteria:
    - assemblies_must_have_x_matches (float): Minimum number of distinct query protein matches for a genome to be considered. <1 is a fraction of the number of query proteins, >1 is the number of query proteins.
    - nucleotide_sequences_must_have_x_matches (float): Minimum number of distinct query protein matches required for a nucleotide sequence to be considered. <1 is a fraction of the number of query proteins, >1 is the number of query proteins.
9. All remaining loci are then assigned a cluster, simply moving from start to end and "breaking" into a new cluster if the distance between two loci is > `break_bgc_on_gap_of`
10. Resulting clusters are then filtered to those that contain >= `gene_clusters_must_have_x_matches`
11. Remaining clusters are pulled into the SocialGene object
    - All loci within each clister and +/- `target_bgc_padding` base pairs are pulled from the database
    - All proteins within each cluster are pulled from the database
    - All HMM models annotating proteins within each cluster are pulled from the database
13. Each target BGC is compared again to the input BGC (via domain similarity)
14. Links are created between the input BGC and the target BGCs
15. "Groups" are created by assigning a target protein to a single input BGC protein based on best-hit
    - This is for assigning groups within clustermap.js for visualization purposes only
16. Results are serialized to clustermap.js JSON format and written to disk
