from socialgene.neo4j.neo4j import GraphDriver
from socialgene.utils.logging import log


def _add_to_neo4j(statement, **kwargs):
    try:
        summary = (
            GraphDriver().driver.execute_query(
                statement,
                **kwargs,
                database_="neo4j",
            )
        ).summary
        if summary.metadata.get("stats"):
            log.info(
                f"{summary.metadata.get('stats').get('properties-set')} properties modified"
            )
        else:
            log.info("No properties modified")
    except Exception as e:
        log.warning(e.message)


def _run_transaction_function(statement, **kwargs):
    try:
        log.debug(statement)
        with GraphDriver() as db:
            a = db.run(
                statement,
                **kwargs,
                database_="neo4j",
            ).values()
        return a

    except Exception as e:
        log.error(e)


def add_taxonomic_name_to_assembly(rank):
    _add_to_neo4j(
        f"""
            MATCH (a1:assembly)-[:IS_TAXON|TAXON_PARENT*1..]->(t1:taxid {{rank:"{rank}"}})
            SET a1.{rank} = t1.name;
            """
    )


def add_protein_descriptions():
    log.info(
        "Adding first description to non-redundant proteins (excluding those containing hypothetical); this can take some time."
    )
    log.info("This can take some time")
    res = _run_transaction_function(
        """
                // set protein description as a random input protein's description (except "hypothetical")
                MATCH (p1:protein)
                CALL {
                    WITH p1
                    MATCH (:nucleotide)-[e1:ENCODES]->(p1)
                    WHERE e1.description IS NOT NULL AND NOT toLower(e1.description) CONTAINS "hypothetical"
                    RETURN e1 LIMIT 1
                } IN TRANSACTIONS OF 100 ROWS
                SET p1.description = e1.description
                RETURN count(p1)
            """
    )
    log.info(f"Modified: {res} properties")


def add_antismash_regions_as_edges():
    log.info(
        "Parsing 'import/antismash_results.jsonl' as (:nucleotide)-[e1:ENCODES]->(:proteins) e1.properties"
    )
    _run_transaction_function(
        """
                CALL apoc.load.json("import/antismash_results.jsonl")
                    YIELD value AS jsonl
                CALL {
                    with jsonl
                    UNWIND jsonl.records as record
                    UNWIND keys(record) as nucleotide_id
                    WITH record, nucleotide_id, range(0, size(record[nucleotide_id])-1) as cluster_index
                    UNWIND cluster_index as cluster_index2
                    WITH record, record[nucleotide_id][cluster_index2] as cluster, nucleotide_id,  cluster_index2
                    MATCH (n1:nucleotide {external_id:nucleotide_id})-[c1:ENCODES]->()
                    WHERE cluster.start < (c1.start + 1 ) < cluster.end
                    SET c1.antismash_products = record[nucleotide_id][cluster_index2].products
                    SET c1.antismash_region = cluster_index2
                } IN TRANSACTIONS OF 1000 ROWS;
            """
    )


def add_antismash_regions_as_nodes():
    log.info("Adding new constraint for (:gene_cluster {nuid}), if needed")
    _run_transaction_function(
        """
                CREATE CONSTRAINT gene_cluster IF NOT EXISTS
                FOR (n:gene_cluster)
                REQUIRE (n.uid) IS UNIQUE;
        """
    )
    log.info("Parsing 'import/antismash_results.jsonl' as new antismash nodes")
    _run_transaction_function(
        """
        CALL apoc.load.json("import/antismash_results.jsonl")
        YIELD value AS jsonl
        UNWIND jsonl.records as record
        UNWIND keys(record) as nucleotide_id
        WITH record, nucleotide_id, range(0, size(record[nucleotide_id])-1) as cluster_indices
        UNWIND cluster_indices as cluster_index
        WITH record, record[nucleotide_id][cluster_index] as cluster, nucleotide_id,  cluster_index
        MERGE (a:antismash {uid: nucleotide_id + "_antismash_" + cluster_index})
        SET a:gene_cluster
        SET a += {start: cluster.start, end: cluster.end}
        WITH a, cluster, nucleotide_id
        MATCH (n:nucleotide {external_id: nucleotide_id})
        MERGE (a)<-[:PREDICTED_BGC]-(n)
        WITH a, cluster
           CALL {
              WITH a, cluster
              UNWIND keys(cluster.protoclusters) AS protocluster_key
              WITH a, cluster, cluster.protoclusters[protocluster_key] AS protocluster, protocluster_key
              MERGE (product:product {uid: protocluster.product})
              SET product:antismash
              MERGE (category:category {uid: protocluster.category})
              SET category:antismash
              MERGE (product)-[:IS_A]->(category)
              MERGE (a)-[r:IS_A {tool: protocluster.tool, start: protocluster.start, end: protocluster.end, core_start: protocluster.core_start, core_end: protocluster.core_end}]->(product)
           }

            """
    )


def culture_collections_as_nodes_rels():
    _run_transaction_function(
        """
                MATCH (a1:assembly) where a1.culture_collection is not null
                WITH a1, split(a1.culture_collection,":")[0] as cc_agency
                MERGE (cc:culture_collection {uid:cc_agency})
                MERGE (a1)-[:FOUND_IN]->(cc);
            """
    )
    _run_transaction_function(
        """
                MATCH (ccs:culture_collection)
                MATCH (a1:assembly)
                WHERE a1.strain starts with ccs.uid and not a1.uid starts with "BGC" and a1.culture_collection is null
                MERGE (a1)-[:FOUND_IN]->(ccs)
            """
    )


def set_mibig_bgc():
    _run_transaction_function(
        """
                MATCH (a1:assembly)<-[:ASSEMBLES_TO]-(n1:nucleotide)
                WHERE a1.uid STARTS WITH "BGC00" AND n1.external_id STARTS WITH "BGC00"
                SET a1:mibig_bgc;
            """
    )


def fix_mibig_taxonomy():
    _add_to_neo4j(
        """
                // Link mibig BGCs to taxa if they aren't
                MATCH (n:assembly)
                where n.uid starts with "BGC00" and not (n)-[:IS_TAXON]->() and n.db_xref is not null
                MATCH (t1:taxid {uid:split(n.db_xref, "taxon:")[1]})
                MERGE (n)-[:IS_TAXON]->(t1);

            """
    )


def antismash_as_separate_nodes():
    log.info("Adding new constrain for (:antismash_bgc {nuid})")
    _run_transaction_function(
        """
                CREATE CONSTRAINT antismash_bgc IF NOT EXISTS
                FOR (n:antismash_bgc)
                REQUIRE (n.uid, n.region) IS UNIQUE;
        """
    )
    log.info("Creating new antismash nodes and linking them to proteins")
    _run_transaction_function(
        """
                MATCH (p1:protein)<-[e1:ENCODES]-(n1:nucleotide)
                WHERE  e1.antismash_region is not null
                WITH n1, e1.antismash_region as ar, min(e1.start) as minstart, max(e1.end) as maxend
                CALL {
                WITH n1, ar, minstart, maxend
                CREATE (bgc:antismash_bgc {uid: n1.uid, region:ar, start:minstart, end:maxend})
                CREATE (n1)-[:PREDICTED_BGC]->(bgc)
                } IN TRANSACTIONS OF 1000 rows;
            """
    )
