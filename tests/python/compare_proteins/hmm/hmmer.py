from socialgene.compare_proteins.hmm.hmmer import (
    _calculate_mod_score_not_named,
    _create_tuple,
    calculate_mod_score,
    CompareDomains,
)
import os
from socialgene.config import env_vars
from socialgene.base.socialgene import SocialGene
from socialgene.hmm.hmmer import HMMER
import pytest
import pandas as pd

FIXTURE_DIR = os.path.dirname(os.path.realpath(__file__))
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)
FIXTURE_DIR = os.path.join(FIXTURE_DIR, "data", "test_genomes")
gbk_path = os.path.join(FIXTURE_DIR, "lagriamide_mibig_bgc0001946.gbk")


FIXTURE_DIR = os.path.dirname(os.path.realpath(__file__))
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)
FIXTURE_DIR = os.path.join(FIXTURE_DIR, "data")
hmm_path = os.path.join(FIXTURE_DIR, "pks.hmm")

# gbk_path = "/home/chase/Documents/github/kwan_lab/socialgene/sgpy/tests/python/data/test_genomes/lagriamide_mibig_bgc0001946.gbk"
# hmm_path = "/home/chase/Documents/github/kwan_lab/socialgene/sgpy/tests/python/data"


def test_calculate_mod_score_not_named():
    sg_object = SocialGene()
    sg_object.parse(gbk_path)
    protein_id_list = list(sg_object.proteins.keys())
    sg_object.annotate_with_hmmscan(
        protein_id_list=protein_id_list, hmm_directory=hmm_path, cpus=1
    )
    p1 = sg_object.proteins["Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5"]
    p2 = sg_object.proteins["iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh"]
    res = _calculate_mod_score_not_named(p1, p2)
    assert res == (
        "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
        "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
        36,
        50,
        0.64,
        0.69,
        0.98,
    )


def test_calculate_mod_score():
    sg_object = SocialGene()
    sg_object.parse(gbk_path)
    protein_id_list = list(sg_object.proteins.keys())
    sg_object.annotate_with_hmmscan(
        protein_id_list=protein_id_list, hmm_directory=hmm_path, cpus=1
    )
    p1 = sg_object.proteins["Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5"]
    p2 = sg_object.proteins["iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh"]
    res = calculate_mod_score(p1, p2)
    assert (
        str(type(res))
        == "<class 'socialgene.compare_proteins.hmmer.protein_comparison_modscore'>"
    )
    assert res._asdict() == {
        "query": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
        "target": "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
        "query_n_domains": 36,
        "target_n_domains": 50,
        "levenshtein": 0.64,
        "jaccard": 0.69,
        "mod_score": 0.98,
    }


def test_CompareDomains_compare_one_to_one():
    sg_object = SocialGene()
    sg_object.parse(gbk_path)
    protein_id_list = list(sg_object.proteins.keys())
    sg_object.annotate_with_hmmscan(
        protein_id_list=protein_id_list, hmm_directory=hmm_path, cpus=1
    )
    a = CompareDomains()
    p1 = sg_object.proteins["Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5"]
    p2 = sg_object.proteins["iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh"]
    res = a.compare_one_to_one(p1, p2)
    assert res._asdict() == {
        "query": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
        "target": "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
        "query_n_domains": 36,
        "target_n_domains": 50,
        "levenshtein": 0.64,
        "jaccard": 0.69,
        "mod_score": 0.98,
    }


def test_CompareDomains_compare_one_to_many():
    sg_object = SocialGene()
    sg_object.parse(gbk_path)
    protein_id_list = list(sg_object.proteins.keys())
    sg_object.annotate_with_hmmscan(
        protein_id_list=protein_id_list, hmm_directory=hmm_path, cpus=1
    )
    a = CompareDomains()
    p1 = sg_object.proteins["Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5"]
    p2 = list(sg_object.proteins.values())
    a.compare_one_to_many(p1, p2)

    expected = pd.DataFrame(
        {
            "query": [
                "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            ],
            "target": [
                "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
                "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
                "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
                "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
                "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
            ],
            "score": [1.50, 0.98, 0.85, 0.54, 0.33, 0.13],
        }
    )
    pd.testing.assert_frame_equal(a.df, expected)


def test_CompareDomains_compare_many_to_many():
    sg_object = SocialGene()
    sg_object.parse(gbk_path)
    protein_id_list = list(sg_object.proteins.keys())
    sg_object.annotate_with_hmmscan(
        protein_id_list=protein_id_list, hmm_directory=hmm_path, cpus=1
    )
    a = CompareDomains()
    p1 = list(sg_object.proteins.values())
    p2 = list(sg_object.proteins.values())
    a.compare_many_to_many(p1, p2)
    expected = pd.DataFrame(
        {
            "query": {
                0: "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
                1: "DTee9G4M8sEfnM4HaPfI37rT74pq7M_G",
                2: "9stFB1fGjCdZVZWHLVI3OD4A_DV3WcV6",
                3: "mB22-i4RqtslyO7_HappM4rJ4Z2Qbkfn",
                4: "Qi23auOUTcBzWTmDHuGinrzIuqH7-zVn",
                5: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                6: "DTee9G4M8sEfnM4HaPfI37rT74pq7M_G",
                7: "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
                8: "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
                9: "Tdc2m3PRLsyEzjwyux6BF4arDy2mQ_Bl",
                10: "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
                11: "Qi23auOUTcBzWTmDHuGinrzIuqH7-zVn",
                12: "4_8182J88axMDpFJBZI6kLNJAu8Ittm3",
                13: "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
                14: "IRqRpDzrGB9UhHJD6AzDq_6Xupj00Nte",
                15: "du1Ncfm5UYiFYgDWD8KW1AQJNHlAcVXL",
                16: "WbViYzQw8y-XfCQMgQXkedGduNMJPa14",
                17: "MSHRSCZfdBJP8vdJdaXfeZrThH_4EUMm",
                18: "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
                19: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                20: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                21: "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
                22: "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
                23: "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
                24: "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
                25: "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
                26: "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
                27: "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
                28: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                29: "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
                30: "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
                31: "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
                32: "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
                33: "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
                34: "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
                35: "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
                36: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                37: "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
                38: "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
                39: "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
                40: "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
                41: "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
                42: "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
                43: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                44: "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
                45: "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
                46: "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
                47: "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
            },
            "target": {
                0: "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
                1: "DTee9G4M8sEfnM4HaPfI37rT74pq7M_G",
                2: "9stFB1fGjCdZVZWHLVI3OD4A_DV3WcV6",
                3: "mB22-i4RqtslyO7_HappM4rJ4Z2Qbkfn",
                4: "DTee9G4M8sEfnM4HaPfI37rT74pq7M_G",
                5: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                6: "Qi23auOUTcBzWTmDHuGinrzIuqH7-zVn",
                7: "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
                8: "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
                9: "Tdc2m3PRLsyEzjwyux6BF4arDy2mQ_Bl",
                10: "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
                11: "Qi23auOUTcBzWTmDHuGinrzIuqH7-zVn",
                12: "4_8182J88axMDpFJBZI6kLNJAu8Ittm3",
                13: "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
                14: "IRqRpDzrGB9UhHJD6AzDq_6Xupj00Nte",
                15: "du1Ncfm5UYiFYgDWD8KW1AQJNHlAcVXL",
                16: "WbViYzQw8y-XfCQMgQXkedGduNMJPa14",
                17: "MSHRSCZfdBJP8vdJdaXfeZrThH_4EUMm",
                18: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                19: "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
                20: "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
                21: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                22: "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
                23: "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
                24: "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
                25: "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
                26: "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
                27: "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
                28: "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
                29: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                30: "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
                31: "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
                32: "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
                33: "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
                34: "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
                35: "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
                36: "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
                37: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                38: "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
                39: "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
                40: "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
                41: "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
                42: "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
                43: "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
                44: "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
                45: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                46: "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
                47: "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
            },
            "score": {
                0: 1.5,
                1: 1.5,
                2: 1.5,
                3: 1.5,
                4: 1.5,
                5: 1.5,
                6: 1.5,
                7: 1.5,
                8: 1.5,
                9: 1.5,
                10: 1.5,
                11: 1.5,
                12: 1.5,
                13: 1.5,
                14: 1.5,
                15: 1.5,
                16: 1.5,
                17: 1.5,
                18: 0.98,
                19: 0.98,
                20: 0.85,
                21: 0.85,
                22: 0.77,
                23: 0.77,
                24: 0.6,
                25: 0.6,
                26: 0.55,
                27: 0.55,
                28: 0.54,
                29: 0.54,
                30: 0.47,
                31: 0.47,
                32: 0.36,
                33: 0.36,
                34: 0.35,
                35: 0.35,
                36: 0.33,
                37: 0.33,
                38: 0.24,
                39: 0.24,
                40: 0.2,
                41: 0.2,
                42: 0.13,
                43: 0.13,
                44: 0.13,
                45: 0.13,
                46: 0.09,
                47: 0.09,
            },
        }
    )
    pd.testing.assert_frame_equal(
        a.df.sort_values(["query", "target", "score"], ignore_index=True),
        expected.sort_values(["query", "target", "score"], ignore_index=True),
    )


def test_CompareDomains_compare_all_to_all_parallel():
    sg_object = SocialGene()
    sg_object.parse(gbk_path)
    protein_id_list = list(sg_object.proteins.keys())
    sg_object.annotate_with_hmmscan(
        protein_id_list=protein_id_list, hmm_directory=hmm_path, cpus=1
    )
    a = CompareDomains()
    p2 = list(sg_object.proteins.values())
    a.compare_all_to_all_parallel(p2)
