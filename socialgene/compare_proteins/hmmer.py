from itertools import combinations, product
from multiprocessing import Pool
from socialgene.compare_proteins.base_class import CompareProteinsBaseClass

from collections import namedtuple

from socialgene.scoring.scoring import mod_score
import pandas as pd


_create_tuple = namedtuple(
    "protein_comparison_modscore",
    (
        "query",
        "target",
        "query_n_domains",
        "target_n_domains",
        "levenshtein",
        "jaccard",
        "mod_score",
    ),
)


def _calculate_mod_score_not_named(p1, p2):
    """Compare domains between two proteins

    Args:
        p1 (Protein): SocialGene Protein Object
        p2 (Protein): SocialGene Protein Object

    Returns:
        tuple:
    """
    return (
        p1.hash_id,
        p2.hash_id,
        *mod_score(
            p1.domain_vector,
            p2.domain_vector,
        ),
    )


def calculate_mod_score(p1, p2):
    """Compare domains between two proteins

    Args:
        p1 (Protein): SocialGene Protein Object
        p2 (Protein): SocialGene Protein Object

    Returns:
        namedtuple:
    """
    return _create_tuple(*_calculate_mod_score_not_named(p1, p2))


class CompareDomains(CompareProteinsBaseClass):
    def __init__(self):
        self.protein_comparisons = set()

    @property
    def mod_score_df(self):
        return (
            pd.DataFrame(self.protein_comparisons)
            .sort_values(by=["mod_score"], ascending=False)
            .reset_index(inplace=False, drop=True)
        )

    @property
    def df(self):
        return self.mod_score_df.filter(["query", "target", "mod_score"]).rename(
            columns={"mod_score": "score"}
        )

    def compare_one_to_one(self, p1, p2):
        return calculate_mod_score(p1, p2)

    def compare_one_to_many(self, p1_obj, p2_obj_list):
        for i in p2_obj_list:
            temp = calculate_mod_score(p1_obj, i)
            if temp.jaccard > 0:
                self.protein_comparisons.add(temp)

    def compare_many_to_many(self, p1_obj_list, p2_obj_list):
        for i1, i2 in product(p1_obj_list, p2_obj_list):
            temp = calculate_mod_score(i1, i2)
            if temp.jaccard > 0:
                self.protein_comparisons.add(temp)

    def compare_all_to_all(self, p1_obj_list):
        for i1, i2 in combinations(p1_obj_list, 2):
            temp = calculate_mod_score(i1, i2)
            if temp.jaccard > 0:
                temp
                self.protein_comparisons.add(temp)

    def compare_all_to_all_parallel(self, p1_obj_list, cpus=1):
        # have to use _calculate_mod_score_not_named because named tuple can't pickle "protein_comparison_modscore"
        with Pool(cpus) as p:
            for i in p.starmap(
                _calculate_mod_score_not_named,
                combinations(p1_obj_list, 2),
            ):
                # i[5] == jaccard
                if i[5] > 0.001:
                    self.protein_comparisons.add(_create_tuple(*i))
