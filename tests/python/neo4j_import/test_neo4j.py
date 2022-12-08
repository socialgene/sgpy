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
        "find_similar_bgc": "cbc16a83b73273d46765d18dd3083ea0",
        "get_species_and_assembly_from_protein_list": "58b223980d468f9af0c7e3c2e1ea8a05",
        "get_class_taxid_from_assembly_id": "1a9fb17df38c63732f231603940d6508",
        "get_genus_taxid_from_assembly_id": "03291ae3b81bbb1a4a04bc3f557b9df5",
        "get_protein_name": "b00db14440e07626323d63ddccf44b27",
        "get_protein_info": "8c1c50fad5e82455f92b501ab37f02c1",
        "get_blastp_matches": "eb68a5918df74aaafa6797f513dd73bf",
        "get_protein_domains": "bfb63bbba63c82acbf8a4c0362ea4345",
        "get_hmms_given_protein_ids": "8d0ae57f8cfa263962eec16095c8c78b",
        "get_protein_from_assembly": "2963034efcbc2d61b09287f20ceb2b81",
        "get_assembly_from_proteins": "22ad0cd6ebf55075e21667916bc53435",
        "get_assembly_and_locus_from_protein": "0da3710fab0495cd7546339fbf056284",
        "count_matched_query_proteins_per_assembly": "193802d5edb5267babf1377a01e7cec5",
        "get_all_protein_ids": "b5e8f9862b1b0a87d61a3fb3cc39a6a9",
        "get_tuples_of_proteins_with_shared_domains": "f2c38577a99b13e6f234dcac0d35860d",
        "get_hmm_hashes_for_protein": "f89ab2b6860886371b79529ac4e7d089",
        "get_bgcs": "ae229e1d798e54dfb73e0b8552546132",
        "get_hmms_for_domain_class": "791c575c0c9b44553fb14138b5389f85",
        "find_clusterjs_as_pandas_identicals": "90816ed8c9bdea798110243d47dc1f58",
        "query_1": "29ae5b86b40ef82bc2a9e2f06b464853",
        "find_similar_bgc3": "5e67def5305612a5fa4e52542833f711",
        "search_a_single_protein": "bb5ee2c74a04b1db23a7024f489e057f",
        "retrieve_protein_locations": "22d620882e0b2615561310df7696b0c2",
        "temp1": "ff8cfc5c197e3b74fbabdeb00c25b4cc",
        "temp2": "91416e09f6f570c7dfcf49e8b99be5b7",
    }
