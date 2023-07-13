from collections import defaultdict
from datetime import datetime
from socialgene.hmm.hmmer import HMMER
from socialgene.base.socialgene import SocialGene
from rich import inspect
from pathlib import Path
from socialgene.neo4j.neo4j import GraphDriver
from socialgene.utils.logging import log
import pandas as pd
import numpy as np
from copy import deepcopy
from rich.progress import (
    Progress,
    MofNCompleteColumn,
    TextColumn,
    BarColumn,
    TextColumn,
    TimeElapsedColumn,
    TextColumn,
    TimeRemainingColumn,
)

import logging

logging.getLogger("neo4j").setLevel(logging.WARNING)
logging.getLogger().setLevel(logging.INFO)

then = datetime.now()


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
    """Rank input proteins by how many ANNOTATES relationships will have to be traversed

    Args:
        sg_object (SocialGene): SocialGene object
        reduce_by (int, optional): Divide the total input proteins by this amount. Defaults to 3.
        max_out (int, optional): Set a max number of input proteins to search. Defaults to 5.

    Returns:
        DataFrame: pd.DataFrame with columns ["query_prot_uid", "total_outdegree"]
    """
    # Prioritize based in the cumulative outdegree of HMM nodes per input protein
    input_protein_domain_df = []
    doms = set()
    # Create a dictionary of protein to domains (turned into a DF)
    # and a set of domains across all inputs
    for k, v in sg_object.proteins.items():
        for domain in v.domains:
            input_protein_domain_df.append(
                {"query_prot_uid": k, "hmm_uid": domain.hmm_id}
            )
            doms.add(domain.hmm_id)
    input_protein_domain_df = pd.DataFrame(input_protein_domain_df)
    # input_protein_domain_df["n_hmms"] = input_protein_domain_df.groupby("query_prot_uid")[
    #     "query_prot_uid"
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
        "query_prot_uid"
    )["outdegree"].transform("sum")
    input_protein_domain_df = input_protein_domain_df.sort_values(
        by=["total_outdegree"], ascending=[True], na_position="first"
    )
    input_protein_domain_df = input_protein_domain_df[
        ["query_prot_uid", "total_outdegree"]
    ].drop_duplicates()
    # Take the top x of input proteins that have the lowest outdegree
    out_len = int(len(input_protein_domain_df) / reduce_by)
    if out_len > max_out:
        out_len = max_out
    return input_protein_domain_df.head(max_out)


def _find_sim_protein(protein):
    """Use HMM annotations to find similar proteins in a Neo4j Database

    Args:
        protein (Protein): SocialGene Protein object

    Returns:
        DataFrame: pd.DataFrame with columns
    """
    # This could be done for the whole BGC at once, but is done protein by protein to be more atomic
    dv = list(set(protein.domain_vector))
    with GraphDriver() as db:
        res = db.run(
            """
    WITH $doms AS input_protein_domains
    MATCH (prot1:protein)<-[a1:ANNOTATES]-(h0:hmm)
    // Search 3 HMMs, if present, compare to all
    WHERE h0.uid IN input_protein_domains[0..3]
    WITH input_protein_domains,prot1, count(DISTINCT(h0)) as initial_count
    // size(input_protein_domains[0..3]) in case less than 3 input domains
    WHERE initial_count = size(input_protein_domains[0..3])
    MATCH (prot1)<-[a2:ANNOTATES]-(h1:hmm)
    WITH input_protein_domains, prot1, collect(DISTINCT(h1.uid)) as target_uids
    WHERE apoc.coll.isEqualCollection(input_protein_domains, target_uids)
    MATCH (n1:nucleotide)-[e1:ENCODES]->(prot1)
    RETURN DISTINCT n1.uid as nuc_uid, prot1.uid as target_prot_uid, e1.start as n_start, e1.end as n_end
             """,
            doms=dv,
        )
        return res.data()


def search_for_similar_proteins(input_protein_domain_df, sg_object):
    """Loop _find_sim_protein() through input proteins and get df of all results

    Args:
        input_protein_domain_df (DataFrame()): datafram with
        sg_object (_type_): _description_

    Returns:
        _type_: _description_
    """
    bb = []
    progress_bar = Progress(
        TextColumn("Searching input proteins..."),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        BarColumn(),
        MofNCompleteColumn(),
        TextColumn("• Time elapsed "),
        TimeElapsedColumn(),
    )
    with progress_bar as pg:
        task = pg.add_task("Progress...", total=len(input_protein_domain_df))
        for prot in list(input_protein_domain_df["query_prot_uid"]):
            bb.extend(
                [
                    {"query_prot_uid": prot} | i
                    for i in _find_sim_protein(sg_object.proteins[prot])
                ]
            )
            pg.update(task, advance=1)
    bb = pd.DataFrame(bb)
    zz = bb.sort_values(by=["n_start"], ascending=[True], na_position="first")
    zz.reset_index(inplace=True, drop=True)
    return zz


if not check_for_hmm_outdegree():
    set_hmm_outdegree()

input_protein_domain_df = prioritize_input_proteins(sg_object, reduce_by=1, max_out=50)

z = search_for_similar_proteins(input_protein_domain_df, sg_object)


prioritized_nucleotide_seqs = (
    z.filter(["nuc_uid", "query_prot_uid"])
    .drop_duplicates()
    .groupby(["nuc_uid"])
    .count()
    .sort_values(by=["query_prot_uid"], ascending=[False], na_position="first")
)


prioritized_nucleotide_seqs.reset_index(inplace=True, drop=False)
prioritized_nucleotide_seqs.rename(
    columns={"nuc_uid": "nucleotide_uid", "query_prot_uid": "count_of_matched_inputs"},
    inplace=True,
)

prioritized_nucleotide_seqs = prioritized_nucleotide_seqs[
    prioritized_nucleotide_seqs["count_of_matched_inputs"]
    > len(input_protein_domain_df["query_prot_uid"]) * 0.5
]

prioritized_nucleotide_seqs


############################################################
############################################################


def prioritize_cluster_windows(df, tolerance=40000):
    """Scan a DataFrame representing a single nucleotide sequence and
    find the set of genes within a nucleotide distance of each other that contain the most
    diverse set of hits to input proteins

    Args:
        df (DataFrame): ['query_prot_uid', 'nuc_uid', 'target_prot_uid', 'n_start', 'n_end']
        tolerance (int, optional): Max nucleotide (base pair) distance two genes can be from each other. (some PKS/NRPS can be 10's of thousands bp long) Defaults to 40000.

    Returns:
        _type_: _description_
    """
    # Within an input input_protein_domain_df, find groups of genes within a tolerance of the number of
    # AAs from each other's starts
    # and return a input_protein_domain_df representing the group with the hits to the largest number of distinct input proteins
    s = 0
    e = 0
    max_match = (0, 0, 0)
    # note: [int(i) for i in (df.n_start - df.n_end.shift(1))[1:]]
    # calculates the lagged difference in start/ends of proteins, which necessitates the table
    # having previously being sorted by start positions
    for i, diff in enumerate([int(i) for i in (df.n_start - df.n_end.shift(1))[1:]]):
        if diff > tolerance:
            # len(set(df.iloc[s:e].query_prot_uid)))  count how many input proteins matched within a window
            if (matches := len(set(df.iloc[s:e].query_prot_uid))) > max_match[2]:
                # e+2 to account for np.diff offset and python's half-open intervals
                max_match = (s, e + 2, matches)
            # start new window
            s = i + 1
            e = s
        else:
            # expand window
            e = i
    # account for no splits
    if max_match == (0, 0, 0):
        max_match = (s, e + 2, len(set(df.iloc[s:e].query_prot_uid)))
    if max_match[2] > 2:
        print(max_match)
        return df.iloc[max_match[0] : max_match[1]]


bro = []
for name, group in z[
    z["nuc_uid"].isin(list(prioritized_nucleotide_seqs.nucleotide_uid))
].groupby(["nuc_uid"]):
    bro.append(find_best_cluster(group))






cm_object = Clustermap()
sg_object = SocialGene()
cm_object.create_clustermap_uuids(sg_object=sg_object)
cm_object.add_cluster(sg_object=sg_object)
cm_object.add_groups(sg_object=sg_object)
cm_object.add_links(sg_object=sg_object)



sg_new = SocialGene()

for i in set(list(bro2.query_prot_uid) + list(bro2.target_prot_uid)):
    sg_new.add_protein(hash_id=i)


_ = sg_new.annotate_proteins_with_neo4j(annotate_all=True)




sg_new.compare_proteins(append=True)


_mod_return(sg_new.proteins["4_8182J88axMDpFJBZI6kLNJAu8Ittm3"],sg_new.proteins["FRw07slqpjn6dAN4sRJ9e-8uz652J1m2"])


for i1, i2 in combinations(self.proteins.items(), r=2):
    _mod_return(i1=i1, i2=i2)


sg_object.query_neo4j_for_related_proteins()



###################
    cm_object = Clustermap()
    cm_object.create_clustermap_uuids(sg_object=aa)
    cm_object.add_cluster(sg_object=aa)
    cm_object.add_groups(sg_object=aa, cutoff=1)
    cm_object.add_links(sg_object=aa)
    cm_object.write_clustermap(outpath)
    print(f"{taxon} done")
