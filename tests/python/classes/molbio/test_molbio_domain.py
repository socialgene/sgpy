from collections import OrderedDict
import pytest
from socialgene.base.molbio import Domain
from socialgene.config import env_vars


base_dict = {
    "hmm_id": "hmm_id",
    "env_from": 1,
    "env_to": 100,
    "seq_pro_score": 1.1,
    "evalue": 1.1,
    "domain_bias": 1.1,
    "domain_score": 1.1,
    "seq_pro_bias": 1.1,
    "hmm_from": 1,
    "hmm_to": 100,
    "ali_from": 1,
    "ali_to": 100,
}


temp = (
    [{"exponentialized": False, "i_evalue": 1000}, False],
    [{"exponentialized": False, "i_evalue": 10}, False],
    [{"exponentialized": False, "i_evalue": 1.1}, False],
    [{"exponentialized": False, "i_evalue": 1}, False],
    [{"exponentialized": False, "i_evalue": 0.1}, False],
    [{"exponentialized": False, "i_evalue": 0}, False],
)

combos = [(base_dict | i[0], i[1]) for i in temp]


@pytest.mark.parametrize("a,expected", combos)
def test_domain_1(a, expected):
    env_vars["HMMSEARCH_IEVALUE"] = 0.2
    temp = Domain(**a)
    assert temp.domain_within_threshold() == expected


def test_domain_dict():
    expected = OrderedDict(
        [
            ("hmm_id", "hmm_id"),
            ("env_from", 1),
            ("env_to", 100),
            ("seq_pro_score", 1.1),
            ("evalue", 1),
            ("i_evalue", 1),
            ("domain_bias", 1.1),
            ("domain_score", 1.1),
            ("seq_pro_bias", 1.1),
            ("hmm_from", 1),
            ("hmm_to", 100),
            ("ali_from", 1),
            ("ali_to", 100),
        ]
    )
    env_vars["HMMSEARCH_IEVALUE"] = 0.2
    temp = base_dict | {"exponentialized": False, "i_evalue": 1.1}
    temp = Domain(**temp)
    assert temp.get_dict() == expected
