import os

import pytest

from socialgene.base.molbio import Protein
from socialgene.base.socialgene import SocialGene
from socialgene.compare_proteins.hmm_scoring import _mod_score_tupler, mod_score

FIXTURE_DIR = os.path.dirname(os.path.realpath(__file__))
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)
FIXTURE_DIR = os.path.join(FIXTURE_DIR, "data", "test_genomes")
gbk_path = os.path.join(FIXTURE_DIR, "lagriamide_mibig_bgc0001946.gbk")


FIXTURE_DIR = os.path.dirname(os.path.realpath(__file__))
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)
FIXTURE_DIR = os.path.join(FIXTURE_DIR, "data", "hmms")
hmm_dir = FIXTURE_DIR


def create_mod_score_tupler():
    return _mod_score_tupler(Protein(uid="a"), Protein(uid="b"), *[i for i in range(5)])


@pytest.mark.parametrize("n", [0, 1, 2, 3, 4, 5, 6])
def test_mod_score_tupler(n):
    with pytest.raises(TypeError):
        _mod_score_tupler(*[i for i in range(n)])


def test_create_tuple_type():
    a = create_mod_score_tupler()
    assert isinstance(a, _mod_score_tupler)


def test_create_tuple_fields():
    a = create_mod_score_tupler()
    assert a.__slots__ == (
        "query",
        "target",
        "query_n_domains",
        "target_n_domains",
        "levenshtein",
        "jaccard",
        "mod_score",
    )


def test_create_tuple_repr():
    a = create_mod_score_tupler()
    assert a.__repr__() == (
        "_mod_score_tupler(query=Protein(uid='a'), target=Protein(uid='b'), query_n_domains=0, target_n_domains=1, levenshtein=2, jaccard=3, mod_score=4)"
    )


def test_create_tuple_dict():
    a = create_mod_score_tupler()
    assert a.__dict__() == {
        "query": "a",
        "target": "b",
        "query_n_domains": 0,
        "target_n_domains": 1,
        "levenshtein": 2,
        "jaccard": 3,
        "mod_score": 4,
    }


def test_mod_score():
    sg_object = SocialGene()
    sg_object.parse(gbk_path)
    protein_id_list = list(sg_object.proteins.keys())
    sg_object.annotate_proteins_with_hmmscan(
        protein_id_list=protein_id_list,
        hmm_filepath=os.path.join(hmm_dir, "pks.hmm"),
        cpus=1,
    )
    sg_object.annotate_proteins_with_hmmscan(
        protein_id_list=protein_id_list,
        hmm_filepath=os.path.join(hmm_dir, "single_hmm.hmm"),
        cpus=1,
    )
    p1 = sg_object.proteins["Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5"]
    p2 = sg_object.proteins["iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh"]
    res = mod_score(p1, p2)
    assert res.jaccard == 1
    assert res.levenshtein == 0.72
    assert res.mod_score == 1.23
    assert res.query_n_domains == 31
    assert res.target_n_domains == 40
    assert res.query == p1
    assert res.target == p2


def test_mod_score_no_domains():
    sg_object = SocialGene()
    sg_object.parse(gbk_path)
    p1 = sg_object.proteins["Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5"]
    p2 = sg_object.proteins["iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh"]
    res = mod_score(p1, p2)
    assert res.jaccard == 0
    assert res.levenshtein == 0
    assert res.mod_score == 0
    assert res.query_n_domains == 0
    assert res.target_n_domains == 0
    assert res.query == p1
    assert res.target == p2


def test_mod_score_same_hash_with_domains():
    sg_object = SocialGene()
    sg_object.parse(gbk_path)
    p1 = sg_object.proteins["Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5"]
    p2 = sg_object.proteins["Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5"]
    protein_id_list = list(sg_object.proteins.keys())
    sg_object.annotate_proteins_with_hmmscan(
        protein_id_list=protein_id_list,
        hmm_filepath=os.path.join(hmm_dir, "pks.hmm"),
        cpus=1,
    )
    sg_object.annotate_proteins_with_hmmscan(
        protein_id_list=protein_id_list,
        hmm_filepath=os.path.join(hmm_dir, "single_hmm.hmm"),
        cpus=1,
    )
    res = mod_score(p1, p2)
    assert res.jaccard == 1
    assert res.levenshtein == 1
    assert res.mod_score == 1.5
    assert res.query_n_domains == 31
    assert res.target_n_domains == 31
    assert res.query == p1
    assert res.target == p2


def test_mod_score_same_hash_no_domains():
    sg_object = SocialGene()
    sg_object.parse(gbk_path)
    p1 = sg_object.proteins["Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5"]
    p2 = sg_object.proteins["Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5"]
    res = mod_score(p1, p2)
    assert res.jaccard == 1
    assert res.levenshtein == 1
    assert res.mod_score == 1.5
    assert res.query_n_domains == 0
    assert res.target_n_domains == 0
    assert res.query == p1
    assert res.target == p2


def test_mod_score_input_error():
    with pytest.raises(TypeError):
        mod_score("p1", 1)
