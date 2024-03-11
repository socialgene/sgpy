from abc import ABC, abstractmethod
from typing import List

import pandas as pd

from socialgene.base.molbio import Protein

BlastTab_COLUMNS = {
    "query": str,
    "target": str,
    "pident": "Float64",
    "length": "Int64",
    "mismatch": "Int64",
    "gapopen": "Int64",
    "qstart": "Int64",
    "qend": "Int64",
    "sstart": "Int64",
    "send": "Int64",
    "evalue": "Float64",
    "score": "Float64",
    "qlen": "Int64",
    "slen": "Int64",
}


class ProteinComparison(ABC):
    @abstractmethod
    def compare_proteins(self, p1: List[Protein], p2: List[Protein]): ...  # noqa: E704

    def reciprocal_hits(
        self,
        q_vs_t_df: pd.DataFrame,
        t_vs_q_df: pd.DataFrame,
    ):
        """This function takes a dataframe and returns a dataframe with the reciprocal best hits"""
        # Filter the DF, keep the best target hit for each query protein
        q_vs_t_df = self.best_hit_to_query(q_vs_t_df)
        # Filter the DF, keep the best query hit for each target protein
        t_vs_q_df = self.best_hit_to_query(t_vs_q_df)
        t_vs_q_df.rename(columns={"query": "target", "target": "query"}, inplace=True)
        df = pd.concat([q_vs_t_df, t_vs_q_df])
        del q_vs_t_df, t_vs_q_df
        df.sort_values(by=self.score_column, ascending=False, inplace=True)
        # remove non reciprocal
        df = df[df.duplicated(["query", "target"], keep=False)]
        df = df.drop_duplicates(subset=["query", "target", self.score_column])
        # keep the reciprocal hit with the highest bitscore
        # but if multiple hits have the same highest bitscore, keep them all
        return df[
            df.groupby(["query", "target"])[self.score_column].transform("max")
            == df[self.score_column]
        ]

    def best_hit_to_query(
        self,
        df: pd.DataFrame,
    ):
        return (
            df.reset_index(inplace=False, drop=True)
            .sort_values(self.score_column, ascending=False)
            .drop_duplicates("query", keep="first")
            .reset_index(inplace=False, drop=True)
        )


class HMMDataFrame(ProteinComparison):
    def __init__(self, score_column="score"):
        self.score_column = score_column


class BlastTab(ProteinComparison):
    def __init__(self, score_column="score"):
        self.score_column = score_column
        self.score_scale = float("inf")
