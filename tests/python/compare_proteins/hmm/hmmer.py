from socialgene.compare_proteins.hmm.hmmer import (
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
        "jaccard": 1,
        "mod_score": 1.14,
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
            "query": {
                0: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                1: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                2: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                3: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                4: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                5: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            },
            "target": {
                0: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                1: "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
                2: "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
                3: "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
                4: "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
                5: "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
            },
            "score": {0: 1.5, 1: 1.14, 2: 0.78, 3: 0.6, 4: 0.52, 5: 0.18},
        }
    )
    pd.testing.assert_frame_equal(a.df, expected)
    a.compare_one_to_many(p1, p2, filter=False)
    expected = pd.DataFrame(
        {
            "query": {
                0: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                1: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                2: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                3: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                4: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                5: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                6: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                7: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                8: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                9: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                10: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                11: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                12: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                13: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                14: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                15: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                16: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                17: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                18: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                19: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                20: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                21: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            },
            "target": {
                0: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                1: "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
                2: "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
                3: "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
                4: "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
                5: "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
                6: "Qi23auOUTcBzWTmDHuGinrzIuqH7-zVn",
                7: "5IYMhENey2WCMrKPUz3AqBIZuFSv6DPP",
                8: "mB22-i4RqtslyO7_HappM4rJ4Z2Qbkfn",
                9: "MSHRSCZfdBJP8vdJdaXfeZrThH_4EUMm",
                10: "WbViYzQw8y-XfCQMgQXkedGduNMJPa14",
                11: "4_8182J88axMDpFJBZI6kLNJAu8Ittm3",
                12: "du1Ncfm5UYiFYgDWD8KW1AQJNHlAcVXL",
                13: "DTee9G4M8sEfnM4HaPfI37rT74pq7M_G",
                14: "-l7xLyFZbiZENPLq_GML8JyTRF1Srawr",
                15: "nk3UJUyLaWr8LochtohJ9L6eugdChZL9",
                16: "Tdc2m3PRLsyEzjwyux6BF4arDy2mQ_Bl",
                17: "T_DzOorDp3ROhRRBtuXP3xyAPorpTVD0",
                18: "IRqRpDzrGB9UhHJD6AzDq_6Xupj00Nte",
                19: "ptq1NGhBcUp3TIEqvAUxnnp4LOKwINvn",
                20: "ewjVum5PbpEJA4rl-BfnCAypKl5HXb7x",
                21: "9stFB1fGjCdZVZWHLVI3OD4A_DV3WcV6",
            },
            "score": {
                0: 1.5,
                1: 1.14,
                2: 0.78,
                3: 0.6,
                4: 0.52,
                5: 0.18,
                6: 0.0,
                7: 0.0,
                8: 0.0,
                9: 0.0,
                10: 0.0,
                11: 0.0,
                12: 0.0,
                13: 0.0,
                14: 0.0,
                15: 0.0,
                16: 0.0,
                17: 0.0,
                18: 0.0,
                19: 0.0,
                20: 0.0,
                21: 0.0,
            },
        }
    )


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
                0: "9stFB1fGjCdZVZWHLVI3OD4A_DV3WcV6",
                1: "mB22-i4RqtslyO7_HappM4rJ4Z2Qbkfn",
                2: "T_DzOorDp3ROhRRBtuXP3xyAPorpTVD0",
                3: "4_8182J88axMDpFJBZI6kLNJAu8Ittm3",
                4: "nk3UJUyLaWr8LochtohJ9L6eugdChZL9",
                5: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                6: "Tdc2m3PRLsyEzjwyux6BF4arDy2mQ_Bl",
                7: "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
                8: "-l7xLyFZbiZENPLq_GML8JyTRF1Srawr",
                9: "ewjVum5PbpEJA4rl-BfnCAypKl5HXb7x",
                10: "Qi23auOUTcBzWTmDHuGinrzIuqH7-zVn",
                11: "ptq1NGhBcUp3TIEqvAUxnnp4LOKwINvn",
                12: "5IYMhENey2WCMrKPUz3AqBIZuFSv6DPP",
                13: "WbViYzQw8y-XfCQMgQXkedGduNMJPa14",
                14: "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
                15: "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
                16: "DTee9G4M8sEfnM4HaPfI37rT74pq7M_G",
                17: "Qi23auOUTcBzWTmDHuGinrzIuqH7-zVn",
                18: "du1Ncfm5UYiFYgDWD8KW1AQJNHlAcVXL",
                19: "IRqRpDzrGB9UhHJD6AzDq_6Xupj00Nte",
                20: "MSHRSCZfdBJP8vdJdaXfeZrThH_4EUMm",
                21: "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
                22: "DTee9G4M8sEfnM4HaPfI37rT74pq7M_G",
                23: "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
                24: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                25: "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
                26: "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
                27: "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
                28: "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
                29: "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
                30: "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
                31: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                32: "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
                33: "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
                34: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                35: "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
                36: "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
                37: "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
                38: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                39: "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
                40: "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
                41: "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
                42: "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
                43: "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
                44: "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
                45: "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
                46: "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
                47: "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
                48: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                49: "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
                50: "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
                51: "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
                52: "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
                53: "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
            },
            "target": {
                0: "9stFB1fGjCdZVZWHLVI3OD4A_DV3WcV6",
                1: "mB22-i4RqtslyO7_HappM4rJ4Z2Qbkfn",
                2: "T_DzOorDp3ROhRRBtuXP3xyAPorpTVD0",
                3: "4_8182J88axMDpFJBZI6kLNJAu8Ittm3",
                4: "nk3UJUyLaWr8LochtohJ9L6eugdChZL9",
                5: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                6: "Tdc2m3PRLsyEzjwyux6BF4arDy2mQ_Bl",
                7: "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
                8: "-l7xLyFZbiZENPLq_GML8JyTRF1Srawr",
                9: "ewjVum5PbpEJA4rl-BfnCAypKl5HXb7x",
                10: "Qi23auOUTcBzWTmDHuGinrzIuqH7-zVn",
                11: "ptq1NGhBcUp3TIEqvAUxnnp4LOKwINvn",
                12: "5IYMhENey2WCMrKPUz3AqBIZuFSv6DPP",
                13: "WbViYzQw8y-XfCQMgQXkedGduNMJPa14",
                14: "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
                15: "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
                16: "Qi23auOUTcBzWTmDHuGinrzIuqH7-zVn",
                17: "DTee9G4M8sEfnM4HaPfI37rT74pq7M_G",
                18: "du1Ncfm5UYiFYgDWD8KW1AQJNHlAcVXL",
                19: "IRqRpDzrGB9UhHJD6AzDq_6Xupj00Nte",
                20: "MSHRSCZfdBJP8vdJdaXfeZrThH_4EUMm",
                21: "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
                22: "DTee9G4M8sEfnM4HaPfI37rT74pq7M_G",
                23: "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
                24: "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
                25: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                26: "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
                27: "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
                28: "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
                29: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                30: "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
                31: "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
                32: "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
                33: "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
                34: "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
                35: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                36: "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
                37: "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
                38: "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
                39: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                40: "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
                41: "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
                42: "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
                43: "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
                44: "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
                45: "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
                46: "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
                47: "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
                48: "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
                49: "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
                50: "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
                51: "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
                52: "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
                53: "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
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
                18: 1.5,
                19: 1.5,
                20: 1.5,
                21: 1.5,
                22: 1.5,
                23: 1.5,
                24: 1.14,
                25: 1.14,
                26: 0.79,
                27: 0.79,
                28: 0.78,
                29: 0.78,
                30: 0.78,
                31: 0.78,
                32: 0.66,
                33: 0.66,
                34: 0.6,
                35: 0.6,
                36: 0.59,
                37: 0.59,
                38: 0.52,
                39: 0.52,
                40: 0.51,
                41: 0.51,
                42: 0.46,
                43: 0.46,
                44: 0.42,
                45: 0.42,
                46: 0.24,
                47: 0.24,
                48: 0.18,
                49: 0.18,
                50: 0.18,
                51: 0.18,
                52: 0.16,
                53: 0.16,
            },
        }
    )
    pd.testing.assert_frame_equal(
        a.df.sort_values(["query", "target", "score"], ignore_index=True),
        expected.sort_values(["query", "target", "score"], ignore_index=True),
    )
    a.compare_many_to_many(p1, p2, filter=False)
    # 22 inputs * 22 inputs = 484
    assert len(a.df) == 484
    pd.testing.assert_frame_equal(
        pd.DataFrame(a.df["score"].value_counts()).reset_index(),
        pd.DataFrame(
            {
                "score": {
                    0: 0.0,
                    1: 1.5,
                    2: 0.78,
                    3: 0.18,
                    4: 1.14,
                    5: 0.79,
                    6: 0.66,
                    7: 0.6,
                    8: 0.59,
                    9: 0.52,
                    10: 0.51,
                    11: 0.46,
                    12: 0.42,
                    13: 0.24,
                    14: 0.16,
                },
                "count": {
                    0: 430,
                    1: 24,
                    2: 4,
                    3: 4,
                    4: 2,
                    5: 2,
                    6: 2,
                    7: 2,
                    8: 2,
                    9: 2,
                    10: 2,
                    11: 2,
                    12: 2,
                    13: 2,
                    14: 2,
                },
            }
        ),
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
