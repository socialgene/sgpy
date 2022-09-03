import pytest
from socialgene.utils.simple_math import find_exp


def test_find_exp():
    assert find_exp(0.01) == int(-2)
    assert find_exp(1e-2) == int(-2)
    assert find_exp(1e-200) == int(-200)
    assert find_exp(1e200) == int(200)


def test_find_exp_fail():
    with pytest.raises(ValueError):
        find_exp(0)
    with pytest.raises(ValueError):
        find_exp(-1)
