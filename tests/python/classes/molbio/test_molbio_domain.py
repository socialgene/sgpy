import pytest
from socialgene.base.molbio import Domain
from socialgene.config import env_vars


temp = Domain(
    hmm_id="hmm_id",
    env_from=1,
    env_to=100,
    seq_pro_score=1.1,
    evalue=1.1,
    i_evalue=1.1,
    domain_bias=1.1,
    domain_score=1.1,
    seq_pro_bias=1.1,
    hmm_from=1,
    hmm_to=100,
    ali_from=1,
    ali_to=100,
)


def test_domain_1():
    assert temp.get_hmm_id() == "hmm_id"


def test_domain_2():
    env_vars["HMMSEARCH_IEVALUE"] = 0
    assert temp.domain_within_threshold() is True
    env_vars["HMMSEARCH_IEVALUE"] = 0.1
    assert temp.domain_within_threshold() is False
    env_vars["HMMSEARCH_IEVALUE"] = 1
    assert temp.domain_within_threshold() is True
    env_vars["HMMSEARCH_IEVALUE"] = 0.0001
    assert temp.domain_within_threshold() is False


def test_domain_2_fail():
    with pytest.raises(ValueError):
        env_vars["HMMSEARCH_IEVALUE"] = -1
        temp.domain_within_threshold()


temp2 = Domain(
    hmm_id="hmm_id",
    env_from=1,
    env_to=100,
    seq_pro_score=1.1,
    evalue=1.1,
    i_evalue=1e-10,
    domain_bias=1.1,
    domain_score=1.1,
    seq_pro_bias=1.1,
    hmm_from=1,
    hmm_to=100,
    ali_from=1,
    ali_to=100,
)


def test_domain_2_2():
    env_vars["HMMSEARCH_IEVALUE"] = 0
    assert temp2.domain_within_threshold() is True
    env_vars["HMMSEARCH_IEVALUE"] = 1e-11
    assert temp2.domain_within_threshold() is False
    env_vars["HMMSEARCH_IEVALUE"] = 1e-9
    assert temp2.domain_within_threshold() is True
    env_vars["HMMSEARCH_IEVALUE"] = 2
    assert temp2.domain_within_threshold() is True


def test_domain_2_2_fail():
    with pytest.raises(ValueError):
        env_vars["HMMSEARCH_IEVALUE"] = -1
        temp2.domain_within_threshold()
