from typing import Dict, List

from socialgene.neo4j.neo4j import Neo4jQuery
from socialgene.utils.logging import log


def find_similar_proteins(sg_prot: List) -> Dict:
    """Find proteins in the neo4j database that have a list of hmms with a jaccard similarity > x with the input protein

    Args:
        sg_prot (List): list of socialgene Protein class objects (eg `[sg_object.proteins['qEcFCOXPHpf_D6FiGPea8DkV_AWWZjMy']]` )

    Returns:
        Dict: {query_prot_id:[target_prot_ids]}
    """
    log.info(
        "Searching the Neo4j database for proteins with similarity to the input's HMMs"
    )
    temp = Neo4jQuery.query_neo4j(
        cypher_name="find_similar_bgc3",
        param=[
            {
                "prot": k,
                "domains": v.domain_vector(only_unique=True),
            }
            for k, v in sg_prot.items()
        ],
    )
    return temp


def add_query_matches_to_result_sg_object(self):
    for matches in self.query_and_match.values():
        for protein_hash in matches:
            self.result_sg_object.add_protein(hash_id=protein_hash)


def calculate_mod_scores(self):
    for k, v in self.query_and_match.items():
        for i in v:
            self._compare_domain_lists(
                protein_id_1=k,
                protein_id_2=i,
                input_list_1=self.input_sg_object.proteins[k].domain_vector,
                input_list_2=self.result_sg_object.proteins[i].domain_vector,
                append=True,
            )


def query_neo4j_for_related_proteins(protein_dict: Dict) -> List:
    """Search the neo4j database for similar BGCs

    Args:
        protein_dict (Dict): SocialGene class' protein dict

    Returns:
        List: [{"query":input_protein_id, "target":["db_protein_id"]}]
    """
    return Neo4jQuery.query_neo4j(
        cypher_name="find_similar_bgc3",
        param=[
            {"prot": k, "domains": v.domain_vector(only_unique=True)}
            for k, v in protein_dict.items()
        ],
    )
