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
    # TODO: force subclasses to implement self.score_column and self.score_scale and self.name

    @abstractmethod
    def compare_proteins(self, p1: List[Protein], p2: List[Protein]):
        ...

    def reciprocal_hits(
        self,
        q_vs_t_df: pd.DataFrame,
        t_vs_q_df: pd.DataFrame,
    ):
        """This function takes a dataframe and returns a dataframe with the reciprocal best hits"""
        q_vs_t_df = self.best_hit_to_query(q_vs_t_df)
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

    # def reciprocal_hits(self, q_vs_t_df: pd.DataFrame, t_vs_q_df: pd.DataFrame):
    #     # Create a new column in both dataframes: normalised bitscore
    #     q_vs_t_df["norm_bitscore"] = q_vs_t_df.bitscore / q_vs_t_df.length
    #     t_vs_q_df["norm_bitscore"] = t_vs_q_df.bitscore / t_vs_q_df.length
    #     # Create query and subject coverage columns in both dataframes
    #     q_vs_t_df["qcov"] = (q_vs_t_df.length / q_vs_t_df.qlen).clip(upper=1)
    #     t_vs_q_df["qcov"] = (t_vs_q_df.length / t_vs_q_df.qlen).clip(upper=1)
    #     q_vs_t_df["scov"] = (q_vs_t_df.length / q_vs_t_df.slen).clip(upper=1)
    #     t_vs_q_df["scov"] = (t_vs_q_df.length / t_vs_q_df.slen).clip(upper=1)
    #     q_vs_t_df = self.best_hit_to_query(q_vs_t_df)
    #     t_vs_q_df = self.best_hit_to_query(t_vs_q_df)
    #     t_vs_q_df.rename(columns={"query": "target", "target": "query"}, inplace=True)
    #     df = pd.concat([q_vs_t_df, t_vs_q_df])
    #     del q_vs_t_df, t_vs_q_df
    #     df.sort_values(by="bitscore", ascending=False, inplace=True)
    #     # remove non reciprocal
    #     df = df[df.duplicated(["query", "target"], keep=False)]
    #     # keep the reciprocal hit with the highest bitscore
    #     # but if multiple hits have the same highest bitscore, keep them all
    #     return df[
    #         df.groupby(["query", "target"])["bitscore"].transform(max) == df["bitscore"]
    #     ]

    # def best_hit_to_query(self, df: pd.DataFrame):
    #     return (
    #         df.reset_index(inplace=False, drop=True)
    #         .sort_values("bitscore", ascending=False)
    #         .drop_duplicates("query", keep="first")
    #         .reset_index(inplace=False, drop=True)
    #     )
