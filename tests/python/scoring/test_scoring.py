from socialgene.scoring.scoring import mod_score
from collections import OrderedDict
import pytest

# TODO: parameterize to cover more


def test_mod_score0():
    assert mod_score(list(range(0, 5)), list(range(0, 5))) == OrderedDict(
        [("l1", 5), ("l2", 5), ("levenshtein", 0), ("jaccard", 1), ("mod_score", 1.5)]
    )


def test_mod_score1():
    assert mod_score(list(range(0, 15)), list(range(0, 5))) == OrderedDict(
        [
            ("l1", 15),
            ("l2", 5),
            ("levenshtein", 10),
            ("jaccard", 0.33),
            ("mod_score", 0.5),
        ]
    )


def test_mod_score1_check_inverse_arg_order():
    # l1 and l2 will be different but scores should be the same
    assert mod_score(list(range(0, 5)), list(range(0, 15))) == OrderedDict(
        [
            ("l1", 5),
            ("l2", 15),
            ("levenshtein", 10),
            ("jaccard", 0.33),
            ("mod_score", 0.5),
        ]
    )


def test_mod_score2():
    assert mod_score(list(range(0, 15)), list(range(15, 30))) == OrderedDict(
        [
            ("l1", 15),
            ("l2", 15),
            ("levenshtein", 15),
            ("jaccard", 0.0),
            ("mod_score", 0),
        ]
    )


def test_mod_score3():
    assert mod_score(list(range(0, 14)), list(range(15, 30))) == OrderedDict(
        [
            ("l1", 14),
            ("l2", 15),
            ("levenshtein", 15),
            ("jaccard", 0.0),
            ("mod_score", 0),
        ]
    )


def test_mod_score4():
    assert mod_score(list(range(0, 5)), [1, 3]) == OrderedDict(
        [("l1", 5), ("l2", 2), ("levenshtein", 3), ("jaccard", 0.4), ("mod_score", 0.6)]
    )


def test_mod_score_empty_1():
    assert mod_score(list(range(0, 5)), []) == OrderedDict(
        [("l1", 5), ("l2", 0), ("levenshtein", 1), ("jaccard", 0.0), ("mod_score", 0)]
    )


def test_mod_score_empty_1_switch():
    assert mod_score([], list(range(0, 5))) == OrderedDict(
        [("l1", 0), ("l2", 5), ("levenshtein", 1), ("jaccard", 0.0), ("mod_score", 0)]
    )


def test_mod_score_empty_both():
    assert mod_score([], []) == OrderedDict(
        [("l1", 0), ("l2", 0), ("levenshtein", 1), ("jaccard", 0.0), ("mod_score", 0)]
    )


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
