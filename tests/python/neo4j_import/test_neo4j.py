from socialgene.neo4j.neo4j import import_queries
import hashlib


def hasher(x):
    return hashlib.md5(str(x).encode("utf-8")).hexdigest()


query_dict = import_queries()


# this doesn't really test anything, more just a backup to flag if the query code was modified
def test_queries_havent_changed():
    assert {k: hasher(v["query"]) for k, v in query_dict.items()} == {
        "minimal_example": "9f82f76631fb7c02ac491a3d1d40efed",
        "total_node_count": "d8b542d2c70e0ce641d7e50758449db6",
        "total_relationship_count": "c52a5173040dd907ea680a4eeb2ebf20",
        "node_label_count": "f353d52ccbc3ed28727b895ef1ac62d0",
        "relationship_label_count": "9bd3d42808fb5f7bca4b8bf90555c4bd",
        "database_parameters": "ed6127588eb2098f9e03db349ea91eb0",
        "search_protein_hash": "0f2395faff081a61a6962d7e8e06f862",
        "get_assembly_ids": "bf492ee93b9554bcbf11697687fb73e4",
        "find_identical_domain_content": "84d402a5e706bdc8317145f15e153d1e",
        "find_similar_bgc": "fc808f5a40a64b98aed9b9c5ee0ba0ae",
        "get_species_and_assembly_from_protein_list": "408fe750a7915e0f2fd03dd9f1fd319a",
        "get_class_taxid_from_assembly_id": "efca8cdc436155f5302485b728bd6fa2",
        "get_genus_taxid_from_assembly_id": "ebb224d5c7afa7c5726818e5f7116a09",
        "get_protein_name": "35c08233306c77f7a09232e1f49d1ca4",
        "get_protein_info": "d403372669009483267f6bedc4bccad3",
        "get_blastp_matches": "375a383816a40135ee4aa6e38f30f37f",
        "get_protein_domains": "8eaef27bbba305807668a753af94f596",
        "get_hmms_given_protein_ids": "c8c3785e508622a30cde9fdc34ed1bdf",
        "get_protein_from_assembly": "54ece1d2c6b4a5e6d21f2b6a6e39693a",
        "get_assembly_from_proteins": "577a22af7beb0a537ab7cddd8417802f",
        "get_assembly_and_locus_from_protein": "eb474d8e785b0ec72dc14f71ea1e097e",
        "count_matched_query_proteins_per_assembly": "33bbbded79aff0b07ac9b240e9a26e52",
        "get_all_protein_ids": "bfcc99b837214bb1343b661d53c03f6a",
        "get_tuples_of_proteins_with_shared_domains": "3723e1270ab3152d8e9bded6df8a28a4",
        "get_hmm_hashes_for_protein": "d1cf56ddfb82273a734c53f590636c9d",
        "get_bgcs": "d600618f1abd353ba88ed16aa693bad0",
        "get_hmms_for_domain_class": "2e0b307b1bc0f14e5401bb9d37bb3ef5",
        "find_clusterjs_as_pandas_identicals": "ebad171129b335f128dc9f1787b2db2c",
        "query_1": "e7a8c86366c804292c40b690a7df4985",
        "find_similar_bgc3": "5eda15ad30e029019e783fbb5c4d342f",
        "search_a_single_protein": "c0deb6d71a95333786b588860b13a59e",
        "retrieve_protein_locations": "3560e25e680b32be6b1b58d2f82da525",
        "temp1": "6493978af72337426fdf692b7069e21f",
        "temp2": "041ddb0984ef53241d9bcc0f50e3d44e",
        "get_mmseqs_cluster_members": "1c1fa9fb33475f2e2e49223009d861b1",
    }
