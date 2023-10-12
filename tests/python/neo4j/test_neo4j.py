import hashlib

from socialgene.neo4j.neo4j import import_queries


def hasher(x):
    return hashlib.md5(str(x).encode("utf-8")).hexdigest()


query_dict = import_queries()


# this doesn't really test anything, more just a backup to flag if the query code was modified
def test_queries_havent_changed():
    assert {k: hasher(v["query"]) for k, v in query_dict.items()} == {
        "minimal_example": "de2904766c2f47b626720e36c228dacd",
        "total_node_count": "d8b542d2c70e0ce641d7e50758449db6",
        "total_relationship_count": "c52a5173040dd907ea680a4eeb2ebf20",
        "node_label_count": "f353d52ccbc3ed28727b895ef1ac62d0",
        "relationship_label_count": "9bd3d42808fb5f7bca4b8bf90555c4bd",
        "database_parameters": "ed6127588eb2098f9e03db349ea91eb0",
        "search_protein_hash": "bd2ceff17ac986772865400d0641945c",
        "get_assembly_ids": "4b5bf5a52a0595576d973ff0d2393d57",
        "find_identical_domain_content": "750f1128bcfb21e8712b76692900fa17",
        "find_similar_bgc": "502ed18c220370296909a240e455202d",
        "get_species_and_assembly_from_protein_list": "2b3e9e53eaa33c6dd94e3653358083a0",
        "get_class_taxid_from_assembly_id": "41c6e4cee9b7c5061aa121fee94f705b",
        "get_genus_taxid_from_assembly_id": "29cba076fab25e15c416690389edab53",
        "get_protein_name": "52ca77cd56297513b4b362983024ff23",
        "get_protein_info": "1e9b273a2df308c81cd117239d95c0b4",
        "get_blastp_matches": "c329422e30256639a6126a3db1c41e9b",
        "get_protein_domains": "d3c92a9adcedfa885359857468a00151",
        "get_hmms_given_protein_ids": "ac2cdbc1964c833c9f587ba6ec46d1cc",
        "get_protein_from_assembly": "a3744c36f492c18c05d6e750937f44a5",
        "get_assembly_from_proteins": "b4874635e2d84e0cf35a5d67b323572a",
        "get_assembly_and_locus_from_protein": "e843f54e891dabb717b6b77eaaf9bbe6",
        "count_matched_query_proteins_per_assembly": "3f19fc7a870daf92bb28654bf7f0536c",
        "get_all_protein_ids": "615ac50a2048d582918d7a6bfede050d",
        "get_tuples_of_proteins_with_shared_domains": "534edb7d4ed5738187bfdb5a07ea9946",
        "get_hmm_hashes_for_protein": "071ac011e48ca146d413136a3dbd6650",
        "get_bgcs": "fc1be527577f137c519c87a1f2dd8a52",
        "get_hmms_for_domain_class": "b10aabecc58db2b59c95544b06c9df8e",
        "find_clusterjs_as_pandas_identicals": "77912f9af1f06617c17f7f5319f0cfbd",
        "query_1": "4958119679fcf31d0c5d118690250f8c",
        "find_similar_bgc3": "6b25e434b1e5a9160db3da7cbf40019e",
        "search_a_single_protein": "5fdd7485616b68dd3a239e4412daa7f1",
        "retrieve_protein_locations": "8811ee509eebffc9348fae79f97d9a5e",
        "temp1": "0aae42c7e4934470fd3090705089e130",
        "temp2": "999621f653a7383817de601ea71af706",
        "get_mmseqs_cluster_members": "2f1568bf7e09eba432ffd1a5531d5623",
    }
