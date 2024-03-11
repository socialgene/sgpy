import argparse
import multiprocessing
from functools import partial
from itertools import batched

import pandas as pd
from rdkit import DataStructs

from socialgene.addons.chemistry.nr import ChemicalCompoundNode, TanimotoSimilarity
from socialgene.base.chem import ChemicalCompound
from socialgene.neo4j.neo4j import GraphDriver

cmpd_label = ChemicalCompoundNode.neo4j_label[0]


def create_arg_parser():
    parser = argparse.ArgumentParser(
        description="Calculate chemical similarity between all the chemical nodes in a SocialGene Neo4j Database"
    )
    parser.add_argument(
        "--cpus",
        help="Number of cpus to use, memory usage will increase linearly. Default is 1.",
        required=False,
        default=1,
        type=int,
    )
    parser.add_argument(
        "--threshold",
        help="Tanimoto similarity threshold, values equal to or above this will be included in the database. Default is 0.5.",
        required=False,
        default=0.5,
        type=float,
    )
    return parser


def get_db_inchis():
    with GraphDriver() as db:
        res = db.run(
            f"""
            MATCH (c1:{cmpd_label})
            RETURN c1.inchi as chem
            """,
        ).value()
    return res


def inchi_list_to_compound_dict(x):
    res = {}
    for i in x:
        temp = ChemicalCompound(i)
        res[i] = (temp, temp.node)
    return res


def create_inchi_list_to_compound_dict(inchlist, cpus=1):
    group_size = len(inchlist) // cpus
    b = batched(inchlist, group_size)
    combined_dict = {}
    with multiprocessing.Pool(processes=cpus) as pool:
        for res in pool.imap_unordered(inchi_list_to_compound_dict, b):
            combined_dict.update(res)
    return combined_dict


def calculate_similarity(tupgroup, threshold=0.5):
    chemsim_set = set()
    tup = tupgroup[2]
    for iter in range(tupgroup[0], tupgroup[1]):
        for i, x in enumerate(DataStructs.BulkTanimotoSimilarity(tup[iter], tup)):
            if x >= threshold:
                if tup[iter] != i:
                    score = int(x * 100)
                    temp = (iter, i, score)
                    chemsim_set.add(temp)
    return chemsim_set


def calculate_similarity_parallel(combined_dict, threshold=0.5, cpus=1):
    id_list = list(combined_dict.keys())
    big_vec = [i[0].morgan for i in combined_dict.values()]
    group_size = len(id_list) // cpus
    groups = [(i, i + group_size, big_vec) for i in range(0, len(id_list), group_size)]
    groups[-1] = (groups[-1][0], len(id_list) - 1, big_vec)
    chemsim = set()
    with multiprocessing.Pool(processes=cpus) as pool:
        calculate_similarity_partial = partial(
            calculate_similarity, threshold=threshold
        )
        for i in pool.imap_unordered(calculate_similarity_partial, groups):
            chemsim.update(i)
    chemsim2 = set()
    for i in chemsim:
        # sorted prevents identical forward + reverse relationships
        chemsim2.add(tuple(sorted(i[0:2]) + [i[2]]))
    return chemsim2


def combine(combined_dict, chemsim):
    df0 = pd.DataFrame(
        ((k, v[0], v[1]) for k, v in combined_dict.items()),
        columns=["id", "compound", "node"],
    )
    df = pd.DataFrame(chemsim, columns=["start", "end", "similarity"])
    df = df.set_index("start")
    df = df.join(df0["node"])
    df.columns = ["end", "similarity", "start_node"]
    df = df.reset_index(level=0)
    df = df.set_index("end")
    df = df.join(df0["node"])
    df.columns = ["start", "similarity", "start_node", "end_node"]
    df = df.reset_index(level=0)
    df["relationship"] = df.apply(
        lambda x: TanimotoSimilarity(
            start=x.start_node, end=x.end_node, properties={"similarity": x.similarity}
        ),
        axis=1,
    )
    return df


def main():
    parser = create_arg_parser()
    args = parser.parse_args()
    if args.cpus == 1:
        print("Warning: Running with 1 cpu will be slow")
    inchis = get_db_inchis()
    combined_dict = create_inchi_list_to_compound_dict(inchis, cpus=args.cpus)
    chemsim = calculate_similarity_parallel(
        combined_dict, threshold=args.threshold, cpus=args.cpus
    )
    df = combine(combined_dict, chemsim)
    del combined_dict, chemsim
    # remove self relationships
    df = df[df.start != df.end]
    z = df["relationship"].tolist()
    z[0].add_multiple_to_neo4j(z, create=True)


if __name__ == "__main__":
    main()
