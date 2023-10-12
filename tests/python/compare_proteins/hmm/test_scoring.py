import os

import pytest

from socialgene.base.socialgene import SocialGene
from socialgene.compare_proteins.hmm.scoring import _mod_score_tupler, mod_score

FIXTURE_DIR = os.path.dirname(os.path.realpath(__file__))
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)
FIXTURE_DIR = os.path.join(FIXTURE_DIR, "data", "test_genomes")
gbk_path = os.path.join(FIXTURE_DIR, "lagriamide_mibig_bgc0001946.gbk")


FIXTURE_DIR = os.path.dirname(os.path.realpath(__file__))
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)
FIXTURE_DIR = os.path.join(FIXTURE_DIR, "data", "hmms")
hmm_path = FIXTURE_DIR


@pytest.mark.parametrize("n", [0, 1, 2, 3, 4, 5, 6])
def test_mod_score_tupler(n):
    with pytest.raises(TypeError):
        _mod_score_tupler(*[i for i in range(n)])


def test_create_tuple_len():
    assert _mod_score_tupler(*[i for i in range(7)])


def test_create_tuple_type():
    a = _mod_score_tupler(*[i for i in range(7)])
    assert (
        str(type(a))
        == "<class 'socialgene.compare_proteins.hmm.scoring.protein_comparison_modscore'>"
    )


def test_create_tuple_fields():
    a = _mod_score_tupler(*[i for i in range(7)])
    assert a._fields == (
        "query",
        "target",
        "query_n_domains",
        "target_n_domains",
        "levenshtein",
        "jaccard",
        "mod_score",
    )


def test_create_tuple_repr():
    a = _mod_score_tupler(*[i for i in range(7)])
    assert a.__repr__() == (
        "protein_comparison_modscore(query=0, target=1, query_n_domains=2, target_n_domains=3, levenshtein=4, jaccard=5, mod_score=6)"
    )


def test_mod_score():
    sg_object = SocialGene()
    sg_object.parse(gbk_path)
    protein_id_list = list(sg_object.proteins.keys())
    sg_object.annotate_proteins_with_hmmscan(
        protein_id_list=protein_id_list, hmm_directory=hmm_path, cpus=1
    )
    p1 = sg_object.proteins["Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5"]
    p2 = sg_object.proteins["iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh"]
    res = mod_score(p1, p2)
    assert res._asdict() == {
        "query": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
        "target": "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
        "query_n_domains": 36,
        "target_n_domains": 50,
        "levenshtein": 0.64,
        "jaccard": 1,
        "mod_score": 1.14,
    }


def test_mod_score_no_domains():
    sg_object = SocialGene()
    sg_object.parse(gbk_path)
    p1 = sg_object.proteins["Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5"]
    p2 = sg_object.proteins["iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh"]
    res = mod_score(p1, p2)
    assert res._asdict() == {
        "query": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
        "target": "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
        "query_n_domains": 0,
        "target_n_domains": 0,
        "levenshtein": 100,
        "jaccard": 0,
        "mod_score": 0,
    }


def test_mod_score_same_hash_with_domains():
    sg_object = SocialGene()
    sg_object.parse(gbk_path)
    p1 = sg_object.proteins["Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5"]
    p2 = sg_object.proteins["Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5"]
    protein_id_list = list(sg_object.proteins.keys())
    sg_object.annotate_proteins_with_hmmscan(
        protein_id_list=protein_id_list, hmm_directory=hmm_path, cpus=1
    )
    res = mod_score(p1, p2)
    assert res._asdict() == {
        "query": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
        "target": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
        "query_n_domains": 36,
        "target_n_domains": 36,
        "levenshtein": 0,
        "jaccard": 1,
        "mod_score": 1.5,
    }


def test_mod_score_same_hash_no_domains():
    sg_object = SocialGene()
    sg_object.parse(gbk_path)
    p1 = sg_object.proteins["Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5"]
    p2 = sg_object.proteins["Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5"]
    res = mod_score(p1, p2)
    assert res._asdict() == {
        "query": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
        "target": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
        "query_n_domains": 0,
        "target_n_domains": 0,
        "levenshtein": 0,
        "jaccard": 1,
        "mod_score": 1.5,
    }


def test_mod_score_input_error():
    with pytest.raises(TypeError):
        mod_score("p1", 1)
