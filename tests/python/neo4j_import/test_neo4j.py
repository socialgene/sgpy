from socialgene.neo4j.neo4j import import_queries
import hashlib


def hasher(x):
    return hashlib.md5(str(x).encode("utf-8")).hexdigest()


query_dict = import_queries()


# this doesn't really test anything, more just a backup to flag if the query code was modified
def test_queries_havent_changed():
    assert {k: hasher(v["query"]) for k, v in query_dict.items()} == {
        "minimal_example": "8fda6dbe49a1b6a19a2cc530be815e02",
        "total_node_count": "1498bbea6b88e57d5e0dd27c16e0f7c8",
        "total_relationship_count": "c3fe9e1e0597080704542d0a892984f5",
        "node_label_count": "85848555d28d250b808c637e1887cf9f",
        "relationship_label_count": "72e5361c0989257227a0e217a1838d6f",
        "database_parameters": "bf0fc6bbfcb2d52b41224fdd9d73d0f0",
        "search_protein_hash": "f666aaf65441b760563a2e84588c92c9",
        "get_assembly_ids": "5cd2cc837ed334341e74b113b91ead72",
        "find_identical_domain_content": "8818ba0221468da7348869794e915f68",
        "find_similar_bgc": "1f7ea8ac7841a13013c7461a99647507",
        "get_species_and_assembly_from_protein_list": "2599abceffe8bf0978d9e3af53dce79f",
        "get_class_taxid_from_assembly_id": "1a9fb17df38c63732f231603940d6508",
        "get_genus_taxid_from_assembly_id": "03291ae3b81bbb1a4a04bc3f557b9df5",
        "get_protein_name": "b00db14440e07626323d63ddccf44b27",
        "get_protein_info": "8c1c50fad5e82455f92b501ab37f02c1",
        "get_blastp_matches": "eb68a5918df74aaafa6797f513dd73bf",
        "get_protein_domains": "bfb63bbba63c82acbf8a4c0362ea4345",
        "get_hmms_given_protein_ids": "8d0ae57f8cfa263962eec16095c8c78b",
        "get_protein_from_assembly": "a24b73366de142375dff38b8484667cf",
        "get_assembly_from_proteins": "7f61ca742e6021641c8696e55cf84663",
        "get_assembly_and_locus_from_protein": "89c960ad47a98c2a425b3203509a1c81",
        "count_matched_query_proteins_per_assembly": "377e389b0eddd8cd97b9046f961bd464",
        "get_all_protein_ids": "b5e8f9862b1b0a87d61a3fb3cc39a6a9",
        "get_tuples_of_proteins_with_shared_domains": "f2c38577a99b13e6f234dcac0d35860d",
        "get_hmm_hashes_for_protein": "f89ab2b6860886371b79529ac4e7d089",
        "get_bgcs": "792fcd74e16e35777318f57c29755a47",
        "get_hmms_for_domain_class": "791c575c0c9b44553fb14138b5389f85",
        "find_clusterjs_as_pandas_identicals": "1f0003b413d744b4d72f535c298fab69",
        "query_1": "d49f4b9f4d37a7a9c90995c602321263",
        "find_similar_bgc3": "5e67def5305612a5fa4e52542833f711",
        "search_a_single_protein": "bb5ee2c74a04b1db23a7024f489e057f",
        "retrieve_protein_locations": "c2af5849e413b60b0beef045feb5f2a4",
        "temp1": "ff8cfc5c197e3b74fbabdeb00c25b4cc",
        "temp2": "91416e09f6f570c7dfcf49e8b99be5b7",
    }
