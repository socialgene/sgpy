from itertools import combinations, product
from multiprocessing import Pool

import pandas as pd

from socialgene.compare_proteins.base_class import CompareProteinsBaseClass
from socialgene.compare_proteins.hmm.scoring import _mod_score_tupler, mod_score


def picklable_modscore(p1, p2):
    # named tuple doesn't work in multiprocessing
    return mod_score(p1, p2)._asdict()


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
        return mod_score(p1, p2)

    def compare_one_to_many(self, p1_obj, p2_obj_list, filter=True):
        for i in p2_obj_list:
            temp = mod_score(p1_obj, i)
            if not filter or temp.jaccard > 0:
                self.protein_comparisons.add(temp)

    def compare_many_to_many(self, p1_obj_list, p2_obj_list, filter=True):
        for i1, i2 in product(p1_obj_list, p2_obj_list):
            temp = mod_score(i1, i2)
            if not filter or temp.jaccard > 0:
                self.protein_comparisons.add(temp)

    def compare_all_to_all(self, p1_obj_list, filter=True):
        for i1, i2 in combinations(p1_obj_list, 2):
            temp = mod_score(i1, i2)
            if not filter or temp.jaccard > 0:
                self.protein_comparisons.add(temp)

    def compare_all_to_all_parallel(self, p1_obj_list, cpus=1, only_hits=True):
        # have to use _calculate_mod_score_not_named because named tuple can't pickle "protein_comparison_modscore"
        with Pool(cpus) as p:
            for i in p.starmap(
                picklable_modscore,
                combinations(p1_obj_list, 2),
            ):
                if not only_hits or i["jaccard"] > 0.001:
                    self.protein_comparisons.add(_mod_score_tupler(**i))
