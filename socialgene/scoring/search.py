from collections import OrderedDict
from itertools import product
from typing import Dict, List
import asyncio
from neo4j import AsyncGraphDatabase, READ_ACCESS
from socialgene.config import env_vars

import pandas as pd
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    TextColumn,
    TimeElapsedColumn,
)

from socialgene.neo4j.neo4j import GraphDriver
from socialgene.utils.logging import log

progress_bar = Progress(
    TextColumn("Pulling target clusters from database..."),
    TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
    BarColumn(),
    MofNCompleteColumn(),
    TextColumn("• Time elapsed "),
    TimeElapsedColumn(),
)


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
    """
    Writes outdegree onto hmm nodes (not just [:ANNOTATES])
    """
    with GraphDriver() as db:
        _ = db.run(
            """
            MATCH (h1: hmm)
            SET h1.outdegree = apoc.node.degree.out(h1, "ANNOTATES")
            """
        )


def get_outdegree_per_hmm_per_protein(sg_object):
    """
    Get the outdegree for each domain in each protein in an input sg_object
    Args:
        sg_object (SocialGene): SocialGene object
    Returns:
        DataFrame: pd.DataFrame with columns ["protein_uid", "hmm_uid", "outdegree"]
    """
    # Get teh non-redundant set of all domains from all proteins
    doms = set()
    for v in sg_object.proteins.values():
        for i in v.domain_vector:
            doms.add(i)
    if doms == set():
        raise ValueError("No domains in sg_object.proteins")
    doms = list(doms)
    # Retrieve the outdegree for the non-redundant set of domains
    with GraphDriver() as db:
        domain_outdegree_df = db.run(
            """
           WITH $doms as doms
           MATCH (h1:hmm)
           WHERE h1.uid in doms
           RETURN h1.uid as hmm_uid, h1.outdegree as outdegree
         """,
            doms=doms,
        ).to_df()
    domain_outdegree_df["hmm_uid"] = domain_outdegree_df["hmm_uid"].astype("category")
    # Create a DF of proteins and their domains
    temp = pd.DataFrame(
        data=(
            {"protein_uid": k, "hmm_uid": i}
            for k, v in sg_object.proteins.items()
            for i in v.domain_vector
        ),
        dtype="category",
    )
    # return a merged df with columns ["protein_uid", "hmm_uid", "outdegree"]
    return temp.merge(domain_outdegree_df, how="left")


async def _find_simimlar_protein(domain_list, frac: float = 0.75):
    """
    The function `_find_sim_protein` is an asynchronous function that queries a Neo4j graph database to
    find similar proteins based on protein domain information.

    Args:
      protein_domain_dict (Dict[List[str]]): The `protein_domain_dict` parameter is a dictionary where
    the keys are protein uids and the values are lists of domain uids.
      frac (float): The `frac` parameter is a float value that represents the fraction of protein
    domains that need to match in order for a protein to be considered similar. By default, it is set to
    0.75, meaning that at least 75% of the protein domains need to match

    Returns:
      The function `_find_sim_protein` returns a Pandas DataFrame containing the results of the query
    executed in the Neo4j database. The DataFrame has columns `assembly_uid`, `nucleotide_uid`,
    `target_prot_uid`, `n_start`, and `n_end`
    """
    # TODO: move async driver to reg driver class module
    async with AsyncGraphDatabase.driver(
        env_vars["NEO4J_URI"],
        auth=(
            env_vars["NEO4J_USER"],
            env_vars["NEO4J_PASSWORD"],
        ),
    ) as driver:
        res = await driver.execute_query(
            """
            WITH $domain_list AS input_protein_domains
            MATCH (prot1:protein)<-[a1:ANNOTATES]-(h0:hmm)
            WHERE h0.uid IN input_protein_domains
            WITH input_protein_domains, prot1, count(DISTINCT(h0)) as initial_count
            WHERE initial_count > size(input_protein_domains) * $frac
            MATCH (n1:nucleotide)-[e1:ENCODES]->(prot1)
            WHERE (n1)-[:ASSEMBLES_TO]->(:assembly)
            MATCH (a1:assembly)<-[:ASSEMBLES_TO]-(n1)
            RETURN a1.uid as assembly_uid, n1.uid as nucleotide_uid, prot1.uid as target_prot_uid, e1.start as n_start, e1.end as n_end
            """,
            domain_list=list(domain_list),
            frac=frac,
        )
        res = pd.DataFrame(res.records, columns=res.keys)
        return res


async def _find_similar_protein_multiple(dict_of_domain_lists, frac: float = 0.75):
    # create task group
    # TODO: if webserver in future this could be used to control max time of search
    async with asyncio.TaskGroup() as group:
        # create and issue tasks
        tasks = {
            k: group.create_task(_find_simimlar_protein(domain_list=v, frac=frac))
            for k, v in dict_of_domain_lists.items()
        }
    # wait for all tasks to complete...
    # report all results
    # return tasks
    return pd.concat([v.result().assign(query_prot_uid=k) for k, v in tasks.items()])


def search_for_similar_proteins(protein_domain_dict, frac):
    """Loop _find_sim_protein() through input proteins and get df of all results

    Args:
        protein_domain_dict (dict): {protein_uid:[hmm_uids]}
        sg_object (_type_): _description_

    Returns:
        _type_: _description_
    """
    bb = []
    progress_bar = Progress(
        TextColumn("Searching database for input proteins..."),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        BarColumn(),
        MofNCompleteColumn(),
        TextColumn("• Time elapsed "),
        TimeElapsedColumn(),
    )
    with progress_bar as pg:
        task = pg.add_task("Progress...", total=len(protein_domain_dict))
        for k, v in protein_domain_dict.items():
            try:
                bb.extend(
                    [
                        {"query_prot_uid": k} | i
                        for i in _find_sim_protein(domain_list=v, frac=frac)
                    ]
                )
            except Exception:
                pass
            pg.update(task, advance=1)
    return pd.DataFrame(bb)


############################################################
############################################################


def prioritize_input_proteins(sg_object, reduce_by=1, max_out=5):
    """Rank input proteins by how many (:hmm)-[:ANNOTATES]->(:protein) relationships will have to be traversed

    Args:
        sg_object (SocialGene): SocialGene object
        reduce_by (int, optional): Divide the total input proteins by this amount. Defaults to 1, for no reduction.
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
    # Drop duplicate protein,hmm pairs
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
    del z
    # Calculate the total outdegree per protein
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
    return input_protein_domain_df.head(out_len)


def prioritize_cluster_windows(df, tolerance=50000, return_df_not_tuple=False):
    """Scan a DataFrame representing a single nucleotide sequence and
    find the set of genes within a nucleotide distance of each other that contain the most
    diverse set of hits to input proteins

    Args:
        df (DataFrame): ['query_prot_uid', 'nucleotide_uid', 'target_prot_uid', 'n_start', 'n_end']
        tolerance (int, optional): Max nucleotide (base pair) distance two genes can be from each other. (some PKS/NRPS can be 10's of thousands bp long) Defaults to 40000.
        return_df_not_tuple (bool, optional) If True, returns filtered Dataframe instead of tuple

    Returns:
        _type_: _description_
    """
    # Within an input input_protein_domain_df, find groups of genes within a tolerance of the number of
    # AAs from each other's starts
    # and return a input_protein_domain_df representing the group with the hits to the largest number of distinct input proteins
    s = 0
    e = 1
    max_match = (0, 0, 0)
    # note: [int(i) for i in (df.n_start - df.n_end.shift(1))[1:]]
    # calculates the lagged difference in start/ends of proteins, which necessitates the table
    # having previously being sorted by start positions
    for diff in [int(i) for i in (df.n_start - df.n_end.shift(1))[1:]]:
        if diff > tolerance:
            # len(set(df.iloc[s:e].query_prot_uid)))  count how many input proteins matched within a window
            if (matches := len(set(df.iloc[s:e].query_prot_uid))) > max_match[2]:
                max_match = (s, e - 1, matches)
            # start new window
            s = e
            e = s + 1
        else:
            # expand window
            e += 1
    # account for no splits
    if max_match == (0, 0, 0):
        max_match = (s, e - 1, len(set(df.iloc[s:e].query_prot_uid)))
    # +1 because  non-inclusive
    max_match = (max_match[0], max_match[1] + 1, max_match[2])
    if return_df_not_tuple:
        return df.iloc[max_match[0] : max_match[1]]
    else:
        return (
            df.iloc[0].nucleotide_uid,
            df.iloc[max_match[0]].n_start,
            df.iloc[max_match[1]].n_end,
        )


def split_df_into_nucl(df):
    for name, group in df.groupby(["nucleotide_uid"]):
        yield group


# a=["pxchmC6JnOKe3miY4Sl-q2GtE-Bp_EsG", "mtmXphnHA8rXHCnw9NTLx6kfBgS0EH2_", "hntB6KFAW_FBC88h_ddLjN6NC2KEouPX", "ZmesxiO51UP4sObKnCNi9uKAkxKcim6W", "PYaLvUL3QKOPoFRg_Ad9FOs8i_IPUsZT", "KH7Md4jb1JFf_stuDUOdiJt-CbrBJmvA", "EY4gYJguMq0I87toz0hLXqtV-nRDhG-1", "CxeY9ZvRAg1F1da2SsBK3lwmp3A4GpFj", "BWPyB1IrZBw5NKdchqjqqPCDCRwJCLcq", "7MxxSFzI-jhUm3TtyCwoDOBwiL1luO91", "1rz1Hm8TaUfF6oapA2ON38scOvJhLbLM"]
# df.groupby('nucleotide_uid').agg({'n_start':'min', 'n_end':'max'})


def count_matches_per_nucleotide_sequence(df):
    return (
        df.filter(["nucleotide_uid", "query_prot_uid"])
        .drop_duplicates()
        .groupby(["nucleotide_uid"])
        .count()
        .reset_index(drop=False)
        .rename(
            columns={
                "nucleotide_uid": "nucleotide_uid",
                "query_prot_uid": "locus_contains_n_queries",
            },
        )
    )


def count_matches_per_assembly(df):
    with GraphDriver() as db:
        res = db.run(
            """
            WITH $nucid as nuc_ids
            unwind nuc_ids as nuc_id
            MATCH (n1: nucleotide {uid: nuc_id})-[:ASSEMBLES_TO]->(a1:assembly)
            RETURN n1.uid as nucleotide_uid, a1.uid as assembly_uid
            """,
            nucid=list(df["nucleotide_uid"]),
        ).to_df()
    temp = (
        df.merge(res, how="left")
        .filter(["assembly_uid", "query_prot_uid"])
        .drop_duplicates()
        .groupby(["assembly_uid"])
        .count()
        .sort_values(by=["query_prot_uid"], ascending=[False], na_position="first")
        .reset_index(drop=False)
        .rename(
            columns={
                "assembly_uid": "assembly_uid",
                "query_prot_uid": "assembly_contains_n_queries",
            },
        )
    )
    return res.merge(temp, how="inner").drop_duplicates()


def count_matches(df):
    df2 = count_matches_per_nucleotide_sequence(df)
    df3 = count_matches_per_assembly(df)
    return df2.merge(df3, how="inner")
