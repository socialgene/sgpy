# import pandas as pd
# import pytest
# from socialgene.search.hmmer.functions import prioritize_input_proteins

# empty_df = pd.DataFrame({"protein_uid": [], "hmm_uid": [], "outdegree": []})

# empty_df["protein_uid"] = empty_df["hmm_uid"].astype("object")
# empty_df["hmm_uid"] = empty_df["hmm_uid"].astype("object")
# empty_df["outdegree"] = empty_df["outdegree"].astype("int")


# df = pd.DataFrame(
#     {
#         "protein_uid": [
#             "a",
#             "a",
#             "a",
#             "b",
#             "b",
#             "b",
#             "c",
#             "d",
#             "d",
#         ],
#         "hmm_uid": [
#             "a",
#             "a",
#             "b",
#             "a",
#             "b",
#             "c",
#             "a",
#             "b",
#             "b",
#         ],
#         "outdegree": [1, 1, 2, 1, 2, 3, 1, 2, 2],
#     }
# )


# def test_nothing():
#     pd.testing.assert_frame_equal(prioritize_input_proteins(df), df)


# def test_nothing():
#     pd.testing.assert_frame_equal(prioritize_input_proteins(df), df)


# @pytest.mark.parametrize(
#     "test_input,expected",
#     [
#         (
#             1,
#             pd.DataFrame(
#                 {"protein_uid": {6: "c"}, "hmm_uid": {6: "a"}, "outdegree": {6: 1}}
#             ),
#         ),
#         (
#             2,
#             pd.DataFrame(
#                 {
#                     "protein_uid": {0: "a", 1: "a", 2: "a", 6: "c"},
#                     "hmm_uid": {0: "a", 1: "a", 2: "b", 6: "a"},
#                     "outdegree": {0: 1, 1: 1, 2: 2, 6: 1},
#                 }
#             ),
#         ),
#         (
#             3,
#             pd.DataFrame(
#                 {
#                     "protein_uid": {0: "a", 1: "a", 2: "a", 6: "c", 7: "d", 8: "d"},
#                     "hmm_uid": {0: "a", 1: "a", 2: "b", 6: "a", 7: "b", 8: "b"},
#                     "outdegree": {0: 1, 1: 1, 2: 2, 6: 1, 7: 2, 8: 2},
#                 }
#             ),
#         ),
#     ],
# )
# def test_max_proteins(test_input, expected):
#     pd.testing.assert_frame_equal(
#         prioritize_input_proteins(df, max_proteins=test_input), expected
#     )


# @pytest.mark.parametrize(
#     "test_input,expected",
#     [
#         (0, empty_df),
#         (
#             1,
#             pd.DataFrame(
#                 {
#                     "protein_uid": {0: "a", 1: "a", 3: "b", 6: "c"},
#                     "hmm_uid": {0: "a", 1: "a", 3: "a", 6: "a"},
#                     "outdegree": {0: 1, 1: 1, 3: 1, 6: 1},
#                 }
#             ),
#         ),
#     ],
# )
# def test_max_outdegree(test_input, expected):
#     pd.testing.assert_frame_equal(
#         prioritize_input_proteins(df, max_outdegree=test_input), expected
#     )


# def test_max_outdegree_error():
#     with pytest.raises(ValueError):
#         prioritize_input_proteins(df, max_outdegree=0.5)


# def test_max_proteins_error():
#     with pytest.raises(ValueError):
#         prioritize_input_proteins(df, max_proteins=0.5)


# def test_max_domains_error():
#     with pytest.raises(ValueError):
#         prioritize_input_proteins(df, max_domains=0.5)
