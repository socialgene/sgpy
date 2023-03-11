import pytest
from socialgene.base.molbio import Protein
from socialgene.config import env_vars


def test_protein():
    temp = Protein(
        sequence="ARNDCQEGHILKMFPSTWYVXZJU",
        description="description",
        other_id="other_id",
    )
    assert temp.description == "description"
    assert temp.domains == set()
    assert temp.hash_id == "0hMjYRUCOMiDkJnVKlZ4QVMGhG8mkwdb"
    assert temp.other_id == "other_id"
    assert temp.sequence == "ARNDCQEGHILKMFPSTWYVXZJU"


def test_fail():
    with pytest.raises(ValueError):
        Protein()


def test_filter_domains():
    env_vars["HMMSEARCH_IEVALUE"] = 1000
    temp = Protein(
        sequence="ARNDCQEGHILKMFPSTWYVXZJU",
        description="description",
        other_id="other_id",
    )
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
    tempd = (
        {"exponentialized": False, "i_evalue": 3},
        {"exponentialized": False, "i_evalue": 2},
        {"exponentialized": False, "i_evalue": 1},
        {"exponentialized": False, "i_evalue": 0},
        {"exponentialized": False, "i_evalue": 0.1},
    )
    for i in [base_dict | i for i in tempd]:
        temp.add_domain(**i)
    # unfiltered
    print(temp.domains)
    assert len(temp.domains) == 4
    # filtered
    temp.filter_domains()
    assert len(temp.domains) == 4
    env_vars["HMMSEARCH_IEVALUE"] = 100
    temp.filter_domains()
    assert len(temp.domains) == 3
    env_vars["HMMSEARCH_IEVALUE"] = 10
    temp.filter_domains()
    assert len(temp.domains) == 2
    env_vars["HMMSEARCH_IEVALUE"] = 1
    temp.filter_domains()
    assert len(temp.domains) == 1
    env_vars["HMMSEARCH_IEVALUE"] = 0.1
    temp.filter_domains()
    assert len(temp.domains) == 0
