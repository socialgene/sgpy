from abc import ABC, abstractmethod
from collections import namedtuple


class CompareProteinsBaseClass(ABC):
    _create_tuple = namedtuple(
        "standard_protein_comparison",
        (
            "query",
            "target",
            "score",
        ),
    )

    def __init__(self):
        self.protein_comparisons = []

    @property
    @abstractmethod
    def df(self):
        ...
        # return (
        #     pd.DataFrame(self.protein_comparisons)
        #     .sort_values(by=["mod_score"], ascending=False)
        #     .reset_index(inplace=False, drop=True)
        # )

    @abstractmethod
    def compare_one_to_one(self, p1, p2):
        ...

    @abstractmethod
    def compare_one_to_many(self, p1_obj, p2_obj_list):
        ...

    @abstractmethod
    def compare_many_to_many(self, p1_obj_list, p2_obj_list):
        ...

    @abstractmethod
    def compare_all_to_all(self):
        ...
