from collections import OrderedDict

import pytest

from socialgene.base.molbio import Domain, Protein
from socialgene.config import env_vars

BASE_DICT = {
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
    [{"exponentialized": True, "i_evalue": 1000}, False],
    [{"exponentialized": True, "i_evalue": 10}, True],
    [{"exponentialized": True, "i_evalue": 1.1}, True],
    [{"exponentialized": True, "i_evalue": 1}, True],
    [{"exponentialized": True, "i_evalue": 0.1}, True],
    [{"exponentialized": True, "i_evalue": 0}, True],
    [{"exponentialized": False, "i_evalue": 1000}, False],
    [{"exponentialized": False, "i_evalue": 10}, True],
    [{"exponentialized": False, "i_evalue": 1.1}, True],
    [{"exponentialized": False, "i_evalue": 1}, True],
    [{"exponentialized": False, "i_evalue": 0.1}, True],
    [{"exponentialized": False, "i_evalue": 0}, True],
)


combos = [(BASE_DICT | i[0], i[1]) for i in temp]


@pytest.mark.parametrize("a,expected", combos)
def test_domain_1(a, expected):
    env_vars["HMMSEARCH_IEVALUE"] = 100
    temp = Domain(**a)
    assert temp.domain_within_threshold == expected


def test_domain_dict():
    expected = OrderedDict(
        [
            ("hmm_id", "hmm_id"),
            ("env_from", 1),
            ("env_to", 100),
            ("seq_pro_score", 1.1),
            ("evalue", 1.1),
            ("i_evalue", 1.1),
            ("domain_bias", 1.1),
            ("domain_score", 1.1),
            ("seq_pro_bias", 1.1),
            ("hmm_from", 1),
            ("hmm_to", 100),
            ("ali_from", 1),
            ("ali_to", 100),
            ("exponentialized", False),
        ]
    )
    env_vars["HMMSEARCH_IEVALUE"] = 1.5
    temp = BASE_DICT | {"exponentialized": False, "i_evalue": 1.1}
    temp2 = Domain(**temp)
    assert temp2.__dict__ == expected
    temp["exponentialized"] = 1
    with pytest.raises(ValueError):
        Domain(**temp)


def test_fail_add_domain_with_negative():
    env_vars["HMMSEARCH_IEVALUE"] = 1000
    temp = Protein(
        sequence="ARNDCQEGHILKMFPSTWYVXZJU",
        description="description",
        external_protein_id="external_protein_id",
    )
    BASE_DICT = {
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
    tempd = ({"exponentialized": True, "i_evalue": -1},)
    for i in [BASE_DICT | i for i in tempd]:
        with pytest.raises(ValueError):
            temp.add_domain(**i)


def test_dont_fail_add_domain_with_negative():
    env_vars["HMMSEARCH_IEVALUE"] = 1000
    temp = Protein(
        sequence="ARNDCQEGHILKMFPSTWYVXZJU",
        description="description",
        external_protein_id="external_protein_id",
    )
    BASE_DICT = {
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
    tempd = ({"exponentialized": False, "i_evalue": -1},)
    for i in [BASE_DICT | i for i in tempd]:
        temp.add_domain(**i)


def test_dont_add_redundant_domain():
    env_vars["HMMSEARCH_IEVALUE"] = 1000
    temp = Protein(
        sequence="ARNDCQEGHILKMFPSTWYVXZJU",
        description="description",
        external_protein_id="external_protein_id",
    )
    tempd = (
        {"exponentialized": False, "i_evalue": -1},
        {"exponentialized": False, "i_evalue": -1},
    )
    for i in [BASE_DICT | i for i in tempd]:
        temp.add_domain(**i)
    assert len(temp.domains) == 1
    assert temp.domains.pop().__dict__ == {
        "ali_from": 1,
        "ali_to": 100,
        "domain_bias": 1.1,
        "domain_score": 1.1,
        "env_from": 1,
        "env_to": 100,
        "evalue": 1.1,
        "hmm_from": 1,
        "hmm_id": "hmm_id",
        "hmm_to": 100,
        "i_evalue": -1.0,
        "seq_pro_bias": 1.1,
        "seq_pro_score": 1.1,
        "exponentialized": False,
    }
