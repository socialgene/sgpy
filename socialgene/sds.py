from collections import defaultdict
from socialgene.hmm.hmmer import HMMER
from socialgene.base.socialgene import SocialGene
from rich import inspect
from pathlib import Path
from socialgene.neo4j.neo4j import GraphDriver
from socialgene.utils.logging import log
import pandas as pd
import numpy as np
from copy import deepcopy
from rich.progress import Progress


import logging

logging.getLogger("neo4j").setLevel(logging.WARNING)
logging.getLogger().setLevel(logging.INFO)


def check_for_hmm_outdegree():
    """
    Checks a single HMM model to see if it has an outdegree property
    """
    with GraphDriver() as db:
        res = db.run(
            """
           MATCH (h1:hmm)
           return h1.outdegree limit 1
         """
        ).peek()
    if res.value():
        return True
    else:
        return False


def set_hmm_outdegree():
    with GraphDriver() as db:
        _ = db.run(
            """
           MATCH (h1:hmm)-[r:ANNOTATES]->()
           WITH h1, count(r) as deg
           SET h1.outdegree = deg
         """
        )


def prioritize_input_proteins(sg_object, reduce_by=3, max_out=5):
    """"""
    # Prioritize based in the cumulative outdegree of HMM nodes per input protein
    input_protein_domain_df = []
    doms = set()
    # Create a dictionary of protein to domains (turned into a DF)
    # and a set of domains across all inputs
    for k, v in sg_object.proteins.items():
        for domain in v.domains:
            input_protein_domain_df.append({"protein_uid": k, "hmm_uid": domain.hmm_id})
            doms.add(domain.hmm_id)
    input_protein_domain_df = pd.DataFrame(input_protein_domain_df)
    # input_protein_domain_df["n_hmms"] = input_protein_domain_df.groupby("protein_uid")[
    #     "protein_uid"
    # ].transform("count")
    input_protein_domain_df.drop_duplicates(inplace=True)
    # Get the outdegree of the domain models
    with GraphDriver() as db:
        res = db.run(
            """
           WITH $doms as doms
           MATCH (h1:hmm)
           WHERE h1.uid in doms
           RETURN h1.uid as hmm_uid, h1.outdegree as outdegree
         """,
            doms=list(doms),
        )
        z = res.to_df()
    # Calculate the summed outdegree of HMM models for each protein
    input_protein_domain_df = input_protein_domain_df.merge(
        z, left_on="hmm_uid", right_on="hmm_uid"
    )
    input_protein_domain_df["total_outdegree"] = input_protein_domain_df.groupby(
        "protein_uid"
    )["outdegree"].transform("sum")
    input_protein_domain_df = input_protein_domain_df.sort_values(
        by=["total_outdegree"], ascending=[True], na_position="first"
    )
    input_protein_domain_df = input_protein_domain_df[
        ["protein_uid", "total_outdegree"]
    ].drop_duplicates()
    # Take the top x of input proteins that have the lowest outdegree
    out_len = int(len(input_protein_domain_df) / reduce_by)
    if out_len > max_out:
        out_len = max_out
    return input_protein_domain_df.head(max_out)


# def find_sim_protein(protein=sg_object.proteins["M-XQN-LZbGTFgj0rp2dFnIsF9bQm9chC"]):
#     # This could be done for the whole BGC at once, but is done protein by protein to be more atomic
#     dv = list(set(protein.domain_vector))
#     if dv:
#         if len(dv) > 10:
#             tol = 2
#         else:
#             tol = 1
#         with GraphDriver() as db:
#             res = db.run(
#                 """

#         WITH $doms AS input_protein_domains
#     WITH input_protein_domains, size(input_protein_domains) AS input_len
#     MATCH (prot1:protein)<-[a1:ANNOTATES]-(h0:hmm)
#     MATCH (prot1)<-[a2:ANNOTATES]-(h1:hmm)
#     WHERE h0.uid IN input_protein_domains
#     WITH input_len, prot1, count(DISTINCT (h0)) AS hmm_matches,  count(DISTINCT (h1)) AS all_annotations
#     WHERE abs(input_len-hmm_matches) < $tol AND abs(input_len-all_annotations) < $tol
#     MATCH (n1:nucleotide)-[e1:ENCODES]->(prot1)
#     RETURN DISTINCT n1.uid as n, prot1.uid as p, input_len, hmm_matches, all_annotations, e1.start as start, e1.end as end
#              """,
#                 doms=dv,
#                 tol=tol,
#             )
#             return res.data()


# input_protein_domain_df = prioritize_input_proteins(sg_object)


# def search_for_similar_proteins(input_protein_domain_df, sg_object):
#     bb = []
#     for prot in list(input_protein_domain_df["protein_uid"]):
#         bb.extend([{"a": prot} | i for i in find_sim_protein(sg_object.proteins[prot])])
#     bb = pd.DataFrame(bb)
#     zz = bb.sort_values(by=["start"], ascending=[True], na_position="first")
#     zz.reset_index(inplace=True, drop=True)
#     return zz


# prioritized_nucleotide_seqs = (
#     zz.filter(["n", "a"]).drop_duplicates().groupby(["n"]).count()
# )
# prioritized_nucleotide_seqs = prioritized_nucleotide_seqs.sort_values(
#     by=["a"], ascending=[False], na_position="first"
# )
# prioritized_nucleotide_seqs.reset_index(inplace=True, drop=False)
# prioritized_nucleotide_seqs.rename(
#     columns={"n": "nucleotide_uid", "a": "count_of_matched_inputs"}, inplace=True
# )

# prioritized_nucleotide_seqs = prioritized_nucleotide_seqs[
#     prioritized_nucleotide_seqs["count_of_matched_inputs"]
#     > len(input_protein_domain_df["protein_uid"]) * 0.75
# ]

# prioritized_nucleotide_seqs


# if not check_for_hmm_outdegree():
#     set_hmm_outdegree()
