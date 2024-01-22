from itertools import combinations, product
from multiprocessing import Pool

import pandas as pd

from socialgene.compare_proteins.base import HMMDataFrame
from socialgene.compare_proteins.hmm_scoring import _mod_score_tupler, mod_score


def picklable_modscore(p1, p2):
    # named tuple doesn't work in multiprocessing
    return mod_score(p1, p2)._asdict()


class CompareDomains(HMMDataFrame):
    def __init__(self):
        self.name = "HMMER annotation comparison with SocialGene"
        self.score_column = "score"

    def compare_proteins(self, queries, targets, **kwargs):
        return pd.DataFrame(
            (
                {
                    "query": i.query.uid,
                    "target": i.target.uid,
                    "query_n_domains": i.query_n_domains,
                    "target_n_domains": i.target_n_domains,
                    "jaccard": i.jaccard,
                    "score": i.mod_score,
                }
                for i in self.compare_many_to_many(
                    queries, targets, filter_non_hits=True
                )
            ),
        ).drop_duplicates(subset=["query", "target"])

    def compare_one_to_one(self, p1, p2):
        protein_comparisons = set()
        return self.protein_comparisons_df(protein_comparisons)

    def compare_one_to_many(self, p1_obj, p2_obj_list, filter_non_hits=True):
        for i in p2_obj_list:
            temp = mod_score(p1_obj, i)
            if not filter_non_hits or temp.jaccard > 0:
                yield temp

    def compare_many_to_many(self, p1_obj_list, p2_obj_list, filter_non_hits=True):
        for i1, i2 in product(p1_obj_list, p2_obj_list):
            temp = mod_score(i1, i2)
            if not filter_non_hits or temp.jaccard > 0:
                yield temp

    def compare_all_to_all(self, p1_obj_list, filter_non_hits=True):
        for i1, i2 in combinations(p1_obj_list, 2):
            temp = mod_score(i1, i2)
            if not filter_non_hits or temp.jaccard > 0:
                yield temp

    def compare_all_to_all_parallel(self, p1_obj_list, cpus=1, only_hits=True):
        # have to use _calculate_mod_score_not_named because named tuple can't pickle "protein_comparison_modscore"
        with Pool(cpus) as p:
            for i in p.starmap(
                picklable_modscore,
                combinations(p1_obj_list, 2),
            ):
                if not only_hits or i["jaccard"] > 0.001:
                    temp = _mod_score_tupler(**i)
                    yield temp
