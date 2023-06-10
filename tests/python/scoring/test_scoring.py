import itertools
from collections import OrderedDict

import pytest

from socialgene.scoring.scoring import mod_score

# TODO: parameterize to cover more




a = list(
    itertools.chain(
        *list([itertools.combinations(range(0, 5), i) for i in range(0, 6)])
    )
)
a = [list(i) for i in a]

result_list = [
    [("l1", 0), ("l2", 5), ("levenshtein", 0), ("jaccard", 0.0), ("mod_score", 0)],
    [("l1", 1), ("l2", 5), ("levenshtein", 0.2), ("jaccard", 0.2), ("mod_score", 0.3)],
    [("l1", 1), ("l2", 5), ("levenshtein", 0.2), ("jaccard", 0.2), ("mod_score", 0.3)],
    [("l1", 1), ("l2", 5), ("levenshtein", 0.2), ("jaccard", 0.2), ("mod_score", 0.3)],
    [("l1", 1), ("l2", 5), ("levenshtein", 0.2), ("jaccard", 0.2), ("mod_score", 0.3)],
    [("l1", 1), ("l2", 5), ("levenshtein", 0.2), ("jaccard", 0.2), ("mod_score", 0.3)],
    [("l1", 2), ("l2", 5), ("levenshtein", 0.4), ("jaccard", 0.4), ("mod_score", 0.6)],
    [("l1", 2), ("l2", 5), ("levenshtein", 0.4), ("jaccard", 0.4), ("mod_score", 0.6)],
    [("l1", 2), ("l2", 5), ("levenshtein", 0.4), ("jaccard", 0.4), ("mod_score", 0.6)],
    [("l1", 2), ("l2", 5), ("levenshtein", 0.4), ("jaccard", 0.4), ("mod_score", 0.6)],
    [("l1", 2), ("l2", 5), ("levenshtein", 0.4), ("jaccard", 0.4), ("mod_score", 0.6)],
    [("l1", 2), ("l2", 5), ("levenshtein", 0.4), ("jaccard", 0.4), ("mod_score", 0.6)],
    [("l1", 2), ("l2", 5), ("levenshtein", 0.4), ("jaccard", 0.4), ("mod_score", 0.6)],
    [("l1", 2), ("l2", 5), ("levenshtein", 0.4), ("jaccard", 0.4), ("mod_score", 0.6)],
    [("l1", 2), ("l2", 5), ("levenshtein", 0.4), ("jaccard", 0.4), ("mod_score", 0.6)],
    [("l1", 2), ("l2", 5), ("levenshtein", 0.4), ("jaccard", 0.4), ("mod_score", 0.6)],
    [("l1", 3), ("l2", 5), ("levenshtein", 0.6), ("jaccard", 0.6), ("mod_score", 0.9)],
    [("l1", 3), ("l2", 5), ("levenshtein", 0.6), ("jaccard", 0.6), ("mod_score", 0.9)],
    [("l1", 3), ("l2", 5), ("levenshtein", 0.6), ("jaccard", 0.6), ("mod_score", 0.9)],
    [("l1", 3), ("l2", 5), ("levenshtein", 0.6), ("jaccard", 0.6), ("mod_score", 0.9)],
    [("l1", 3), ("l2", 5), ("levenshtein", 0.6), ("jaccard", 0.6), ("mod_score", 0.9)],
    [("l1", 3), ("l2", 5), ("levenshtein", 0.6), ("jaccard", 0.6), ("mod_score", 0.9)],
    [("l1", 3), ("l2", 5), ("levenshtein", 0.6), ("jaccard", 0.6), ("mod_score", 0.9)],
    [("l1", 3), ("l2", 5), ("levenshtein", 0.6), ("jaccard", 0.6), ("mod_score", 0.9)],
    [("l1", 3), ("l2", 5), ("levenshtein", 0.6), ("jaccard", 0.6), ("mod_score", 0.9)],
    [("l1", 3), ("l2", 5), ("levenshtein", 0.6), ("jaccard", 0.6), ("mod_score", 0.9)],
    [("l1", 4), ("l2", 5), ("levenshtein", 0.8), ("jaccard", 0.8), ("mod_score", 1.2)],
    [("l1", 4), ("l2", 5), ("levenshtein", 0.8), ("jaccard", 0.8), ("mod_score", 1.2)],
    [("l1", 4), ("l2", 5), ("levenshtein", 0.8), ("jaccard", 0.8), ("mod_score", 1.2)],
    [("l1", 4), ("l2", 5), ("levenshtein", 0.8), ("jaccard", 0.8), ("mod_score", 1.2)],
    [("l1", 4), ("l2", 5), ("levenshtein", 0.8), ("jaccard", 0.8), ("mod_score", 1.2)],
    [("l1", 5), ("l2", 5), ("levenshtein", 1.0), ("jaccard", 1), ("mod_score", 1.5)],
]

result_list_reversed = [
    [("l1", 5), ("l2", 0), ("levenshtein", 0), ("jaccard", 0.0), ("mod_score", 0)],
    [("l1", 5), ("l2", 1), ("levenshtein", 0.2), ("jaccard", 0.2), ("mod_score", 0.3)],
    [("l1", 5), ("l2", 1), ("levenshtein", 0.2), ("jaccard", 0.2), ("mod_score", 0.3)],
    [("l1", 5), ("l2", 1), ("levenshtein", 0.2), ("jaccard", 0.2), ("mod_score", 0.3)],
    [("l1", 5), ("l2", 1), ("levenshtein", 0.2), ("jaccard", 0.2), ("mod_score", 0.3)],
    [("l1", 5), ("l2", 1), ("levenshtein", 0.2), ("jaccard", 0.2), ("mod_score", 0.3)],
    [("l1", 5), ("l2", 2), ("levenshtein", 0.4), ("jaccard", 0.4), ("mod_score", 0.6)],
    [("l1", 5), ("l2", 2), ("levenshtein", 0.4), ("jaccard", 0.4), ("mod_score", 0.6)],
    [("l1", 5), ("l2", 2), ("levenshtein", 0.4), ("jaccard", 0.4), ("mod_score", 0.6)],
    [("l1", 5), ("l2", 2), ("levenshtein", 0.4), ("jaccard", 0.4), ("mod_score", 0.6)],
    [("l1", 5), ("l2", 2), ("levenshtein", 0.4), ("jaccard", 0.4), ("mod_score", 0.6)],
    [("l1", 5), ("l2", 2), ("levenshtein", 0.4), ("jaccard", 0.4), ("mod_score", 0.6)],
    [("l1", 5), ("l2", 2), ("levenshtein", 0.4), ("jaccard", 0.4), ("mod_score", 0.6)],
    [("l1", 5), ("l2", 2), ("levenshtein", 0.4), ("jaccard", 0.4), ("mod_score", 0.6)],
    [("l1", 5), ("l2", 2), ("levenshtein", 0.4), ("jaccard", 0.4), ("mod_score", 0.6)],
    [("l1", 5), ("l2", 2), ("levenshtein", 0.4), ("jaccard", 0.4), ("mod_score", 0.6)],
    [("l1", 5), ("l2", 3), ("levenshtein", 0.6), ("jaccard", 0.6), ("mod_score", 0.9)],
    [("l1", 5), ("l2", 3), ("levenshtein", 0.6), ("jaccard", 0.6), ("mod_score", 0.9)],
    [("l1", 5), ("l2", 3), ("levenshtein", 0.6), ("jaccard", 0.6), ("mod_score", 0.9)],
    [("l1", 5), ("l2", 3), ("levenshtein", 0.6), ("jaccard", 0.6), ("mod_score", 0.9)],
    [("l1", 5), ("l2", 3), ("levenshtein", 0.6), ("jaccard", 0.6), ("mod_score", 0.9)],
    [("l1", 5), ("l2", 3), ("levenshtein", 0.6), ("jaccard", 0.6), ("mod_score", 0.9)],
    [("l1", 5), ("l2", 3), ("levenshtein", 0.6), ("jaccard", 0.6), ("mod_score", 0.9)],
    [("l1", 5), ("l2", 3), ("levenshtein", 0.6), ("jaccard", 0.6), ("mod_score", 0.9)],
    [("l1", 5), ("l2", 3), ("levenshtein", 0.6), ("jaccard", 0.6), ("mod_score", 0.9)],
    [("l1", 5), ("l2", 3), ("levenshtein", 0.6), ("jaccard", 0.6), ("mod_score", 0.9)],
    [("l1", 5), ("l2", 4), ("levenshtein", 0.8), ("jaccard", 0.8), ("mod_score", 1.2)],
    [("l1", 5), ("l2", 4), ("levenshtein", 0.8), ("jaccard", 0.8), ("mod_score", 1.2)],
    [("l1", 5), ("l2", 4), ("levenshtein", 0.8), ("jaccard", 0.8), ("mod_score", 1.2)],
    [("l1", 5), ("l2", 4), ("levenshtein", 0.8), ("jaccard", 0.8), ("mod_score", 1.2)],
    [("l1", 5), ("l2", 4), ("levenshtein", 0.8), ("jaccard", 0.8), ("mod_score", 1.2)],
    [("l1", 5), ("l2", 5), ("levenshtein", 1.0), ("jaccard", 1), ("mod_score", 1.5)],
]

testdata = list(zip(a, itertools.cycle([[0, 1, 2, 3, 4]]), result_list))
testdatareversed = list(
    zip(a, itertools.cycle([[0, 1, 2, 3, 4]]), result_list_reversed)
)


@pytest.mark.parametrize("a,b,expected", testdata)
def test_mod(a, b, expected):
    score_result = mod_score(a, b)
    assert score_result == OrderedDict(expected)


@pytest.mark.parametrize("a,b,expected", testdatareversed)
def test_mod_reverse(a, b, expected):
    score_result = mod_score(b, a)
    assert score_result == OrderedDict(expected)


# def test_mod_score0():
#     assert mod_score(list(range(0, 5)), list(range(0, 5))) == OrderedDict(
#         [("l1", 5), ("l2", 5), ("levenshtein", 0), ("jaccard", 1), ("mod_score", 1.5)]
#     )


# def test_mod_score1():
#     assert mod_score(list(range(0, 15)), list(range(0, 5))) == OrderedDict(
#         [
#             ("l1", 15),
#             ("l2", 5),
#             ("levenshtein", 0.33),
#             ("jaccard", 0.33),
#             ("mod_score", 0.5),
#         ]
#     )


# def test_mod_score1_check_inverse_arg_order():
#     # l1 and l2 will be different but scores should be the same
#     assert mod_score(list(range(0, 5)), list(range(0, 15))) == OrderedDict(
#         [
#             ("l1", 5),
#             ("l2", 15),
#             ("levenshtein", 10),
#             ("jaccard", 0.33),
#             ("mod_score", 0.5),
#         ]
#     )


# def test_mod_score2():
#     assert mod_score(list(range(0, 15)), list(range(15, 30))) == OrderedDict(
#         [
#             ("l1", 15),
#             ("l2", 15),
#             ("levenshtein", 15),
#             ("jaccard", 0.0),
#             ("mod_score", 0),
#         ]
#     )


# def test_mod_score3():
#     assert mod_score(list(range(0, 14)), list(range(15, 30))) == OrderedDict(
#         [
#             ("l1", 14),
#             ("l2", 15),
#             ("levenshtein", 15),
#             ("jaccard", 0.0),
#             ("mod_score", 0),
#         ]
#     )


# def test_mod_score4():
#     assert mod_score(list(range(0, 5)), [1, 3]) == OrderedDict(
#         [
#             ("l1", 5),
#             ("l2", 2),
#             ("levenshtein", 0.4),
#             ("jaccard", 0.4),
#             ("mod_score", 0.6),
#         ]
#     )


# def test_mod_score_empty_1():
#     assert mod_score(list(range(0, 5)), []) == OrderedDict(
#         [("l1", 5), ("l2", 0), ("levenshtein", 0), ("jaccard", 0.0), ("mod_score", 0)]
#     )


# def test_mod_score_empty_1_switch():
#     assert mod_score([], list(range(0, 5))) == OrderedDict(
#         [("l1", 0), ("l2", 5), ("levenshtein", 1), ("jaccard", 0.0), ("mod_score", 0)]
#     )


# def test_mod_score_empty_both():
#     assert mod_score([], []) == OrderedDict(
#         [("l1", 0), ("l2", 0), ("levenshtein", 1), ("jaccard", 0.0), ("mod_score", 0)]
#     )


def test_mod_score_fail_0():
    with pytest.raises(Exception) as e_info:
        mod_score("a", "a")
        _ = e_info


def test_mod_score_fail_1():
    with pytest.raises(Exception) as e_info:
        mod_score(1, 1)
        _ = e_info


def test_mod_score_fail_2():
    with pytest.raises(Exception) as e_info:
        mod_score([1], 1)
        _ = e_info
