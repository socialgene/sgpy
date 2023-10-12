# import pandas as pd
# import pytest
# from socialgene.search.base import ProcessSearchResult


# @pytest.fixture
# def search_result_df():
#     return pd.DataFrame(
#         {
#             "assembly_uid": ["A", "A", "B", "B", "C", "C"],
#             "nucleotide_uid": ["a", "a", "b", "b", "c", "c"],
#             "target": ["t1", "t2", "t1", "t2", "t1", "t2"],
#             "n_start": [1, 2, 3, 4, 5, 6],
#             "n_end": [10, 20, 30, 40, 50, 60],
#             "query": ["q1", "q2", "q3", "q4", "q5", "q6"],
#         }
#     )


# def test_filter_assemblies(search_result_df):
#     result = ProcessSearchResult(search_result_df, 2)
#     result.filter_assemblies(2)
#     assert len(result.df) == 4
#     assert set(result.df["assembly_uid"]) == {"A", "B"}


# def test_filter_nucleotides(search_result_df):
#     result = ProcessSearchResult(search_result_df, 2)
#     result.filter_nucleotides(2)
#     assert len(result.df) == 4
#     assert set(result.df["nucleotide_uid"]) == {"a", "b"}


# def test_label_clusters(search_result_df):
#     result = ProcessSearchResult(search_result_df, 2)
#     result._label_clusters(max_gap=5)
#     assert len(result.df) == 6
#     assert set(result.df["cluster"]) == {0, 1, 2}


# def test_calc_intrahits(search_result_df):
#     result = ProcessSearchResult(search_result_df, 2)
#     result._label_clusters(max_gap=5)
#     result._calc_intrahits()
#     assert len(result.df) == 6
#     assert set(result.df["cluster_unique_hits"]) == {1, 2}


# def test_process(search_result_df):
#     result = ProcessSearchResult(search_result_df, 2)
#     result.process(2, 2, 5)
#     assert len(result.df) == 4
#     assert set(result.df["assembly_uid"]) == {"A", "B"}
#     assert set(result.df["nucleotide_uid"]) == {"a", "b"}
#     assert set(result.df["cluster"]) == {0, 1}
#     assert set(result.df["cluster_unique_hits"]) == {1, 2}
