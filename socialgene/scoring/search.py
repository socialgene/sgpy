from datetime import datetime
from socialgene.neo4j.neo4j import GraphDriver
from rich.progress import (
    Progress,
    MofNCompleteColumn,
    TextColumn,
    BarColumn,
    TextColumn,
    TimeElapsedColumn,
    TextColumn,
)
import logging
import pandas as pd

logging.getLogger("neo4j").setLevel(logging.WARNING)
logging.getLogger().setLevel(logging.INFO)
then = datetime.now()


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
            SET h1.outdegree = apoc.node.degree(h1)
            """
        )


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
    MATCH (h0:hmm)
    WHERE h0.uid IN input_protein_domains[0..3]
    MATCH (prot1:protein)<-[a1:ANNOTATES]-(h0)
    // Search 3 HMMs, if present, compare to all
    WITH input_protein_domains,prot1, count(DISTINCT(h0)) as initial_count
    // size(input_protein_domains[0..3]) in case less than 3 input domains
    WHERE initial_count = size(input_protein_domains[0..3])
    MATCH (prot1)<-[a2:ANNOTATES]-(h1:hmm)
    WITH input_protein_domains, prot1, collect(DISTINCT(h1.uid)) as target_uids
    WHERE apoc.coll.isEqualCollection(input_protein_domains, target_uids)
    MATCH (n1:nucleotide)-[e1:ENCODES]->(prot1)
    WHERE (n1)-[:ASSEMBLES_TO]->(:assembly)-[:FOUND_IN]->(:culture_collection)
    RETURN DISTINCT n1.uid as nucleotide_uid, prot1.uid as target_prot_uid, e1.start as n_start, e1.end as n_end
             """,
            doms=list(dv),
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
        TextColumn("Searching database for input proteins..."),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        BarColumn(),
        MofNCompleteColumn(),
        TextColumn("• Time elapsed "),
        TimeElapsedColumn(),
    )
    with progress_bar as pg:
        task = pg.add_task("Progress...", total=len(input_protein_domain_df))
        for prot in list(input_protein_domain_df["query_prot_uid"]):
            try:
                bb.extend(
                [
                    {"query_prot_uid": prot} | i
                    for i in _find_sim_protein(sg_object.proteins[prot])
                ]
                )
            except:
                pass
            pg.update(task, advance=1)
    bb = pd.DataFrame(bb)
    zz = bb.sort_values(by=["n_start"], ascending=[True], na_position="first")
    zz.reset_index(inplace=True, drop=True)
    return zz


############################################################
############################################################


def prioritize_cluster_windows(df, tolerance=40000, return_df_not_tuple=False):
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
    # if max_match[2] > 0:
    if return_df_not_tuple:
        return df.iloc[max_match[0] : max_match[1]]
    else:
        return (
            df.iloc[0].nucleotide_uid,
            df.iloc[max_match[0]].n_start,
            df.iloc[max_match[1]].n_end,
        )


def count_matches_per_nucleotide_sequence(df):
    return (
        df.filter(["nucleotide_uid", "query_prot_uid"])
        .drop_duplicates()
        .groupby(["nucleotide_uid"])
        .count()
        .sort_values(by=["query_prot_uid"], ascending=[False], na_position="first")
        .reset_index(drop=False)
        .rename(
            columns={
                "nucleotide_uid": "nucleotide_uid",
                "query_prot_uid": "count_of_matched_inputs",
            },
        )
    )

    