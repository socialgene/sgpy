import os

from socialgene.base.socialgene import SocialGene
from socialgene.compare_proteins.hmmer import CompareDomains

FIXTURE_DIR = os.path.dirname(os.path.realpath(__file__))
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)
FIXTURE_DIR = os.path.join(FIXTURE_DIR, "data", "test_genomes")
gbk_path = os.path.join(FIXTURE_DIR, "lagriamide_mibig_bgc0001946.gbk")


FIXTURE_DIR = os.path.dirname(os.path.realpath(__file__))
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)
FIXTURE_DIR = os.path.join(FIXTURE_DIR, "data", "hmms")
hmm_path = os.path.join(FIXTURE_DIR, "pks.hmm")


def test_CompareDomains_compare_one_to_one():
    sg_object = SocialGene()
    sg_object.parse(gbk_path)
    protein_id_list = list(sg_object.proteins.keys())
    sg_object.annotate_proteins_with_hmmscan(
        protein_id_list=protein_id_list, hmm_filepath=hmm_path, cpus=1
    )
    a = CompareDomains()
    p1 = sg_object.proteins["Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5"]
    p2 = sg_object.proteins["iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh"]
    res = a.compare_one_to_one(p1, p2)
    assert res.to_dict() == {
        "query": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
        "target": "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
        "query_n_domains": 31,
        "target_n_domains": 40,
        "levenshtein": 0.72,
        "jaccard": 1,
        "mod_score": 1.23,
    }


def test_CompareDomains_compare_one_to_many():
    sg_object = SocialGene()
    sg_object.parse(gbk_path)
    protein_id_list = list(sg_object.proteins.keys())
    sg_object.annotate_proteins_with_hmmscan(
        protein_id_list=protein_id_list, hmm_filepath=hmm_path, cpus=1
    )
    a = CompareDomains()
    p1 = sg_object.proteins["Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5"]
    p2 = list(sg_object.proteins.values())
    gen = a.compare_one_to_many(p1, p2)
    expected = [
        {
            "query": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "target": "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
            "query_n_domains": 31,
            "target_n_domains": 14,
            "levenshtein": 0.35,
            "jaccard": 0.46,
            "mod_score": 0.59,
        },
        {
            "query": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "target": "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
            "query_n_domains": 31,
            "target_n_domains": 27,
            "levenshtein": 0.68,
            "jaccard": 0.67,
            "mod_score": 1.01,
        },
        {
            "query": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "target": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "query_n_domains": 31,
            "target_n_domains": 31,
            "levenshtein": 1,
            "jaccard": 1,
            "mod_score": 1.5,
        },
        {
            "query": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "target": "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
            "query_n_domains": 31,
            "target_n_domains": 7,
            "levenshtein": 0.23,
            "jaccard": 0.6,
            "mod_score": 0.53,
        },
        {
            "query": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "target": "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
            "query_n_domains": 31,
            "target_n_domains": 40,
            "levenshtein": 0.72,
            "jaccard": 1,
            "mod_score": 1.23,
        },
        {
            "query": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "target": "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
            "query_n_domains": 31,
            "target_n_domains": 2,
            "levenshtein": 0.06,
            "jaccard": 0.2,
            "mod_score": 0.16,
        },
    ]
    assert [i.to_dict() for i in gen] == expected
    gen = a.compare_one_to_many(p1, p2, filter_non_hits=False)
    expected = [
        {
            "query": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "target": "Tdc2m3PRLsyEzjwyux6BF4arDy2mQ_Bl",
            "query_n_domains": 31,
            "target_n_domains": 2,
            "levenshtein": 0.0,
            "jaccard": 0.0,
            "mod_score": 0,
        },
        {
            "query": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "target": "ptq1NGhBcUp3TIEqvAUxnnp4LOKwINvn",
            "query_n_domains": 31,
            "target_n_domains": 0,
            "levenshtein": 0,
            "jaccard": 0,
            "mod_score": 0,
        },
        {
            "query": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "target": "-l7xLyFZbiZENPLq_GML8JyTRF1Srawr",
            "query_n_domains": 31,
            "target_n_domains": 0,
            "levenshtein": 0,
            "jaccard": 0,
            "mod_score": 0,
        },
        {
            "query": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "target": "T_DzOorDp3ROhRRBtuXP3xyAPorpTVD0",
            "query_n_domains": 31,
            "target_n_domains": 0,
            "levenshtein": 0,
            "jaccard": 0,
            "mod_score": 0,
        },
        {
            "query": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "target": "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
            "query_n_domains": 31,
            "target_n_domains": 14,
            "levenshtein": 0.35,
            "jaccard": 0.46,
            "mod_score": 0.59,
        },
        {
            "query": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "target": "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
            "query_n_domains": 31,
            "target_n_domains": 27,
            "levenshtein": 0.68,
            "jaccard": 0.67,
            "mod_score": 1.01,
        },
        {
            "query": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "target": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "query_n_domains": 31,
            "target_n_domains": 31,
            "levenshtein": 1,
            "jaccard": 1,
            "mod_score": 1.5,
        },
        {
            "query": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "target": "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
            "query_n_domains": 31,
            "target_n_domains": 7,
            "levenshtein": 0.23,
            "jaccard": 0.6,
            "mod_score": 0.53,
        },
        {
            "query": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "target": "DTee9G4M8sEfnM4HaPfI37rT74pq7M_G",
            "query_n_domains": 31,
            "target_n_domains": 1,
            "levenshtein": 0.0,
            "jaccard": 0.0,
            "mod_score": 0,
        },
        {
            "query": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "target": "mB22-i4RqtslyO7_HappM4rJ4Z2Qbkfn",
            "query_n_domains": 31,
            "target_n_domains": 1,
            "levenshtein": 0.0,
            "jaccard": 0.0,
            "mod_score": 0,
        },
        {
            "query": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "target": "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
            "query_n_domains": 31,
            "target_n_domains": 40,
            "levenshtein": 0.72,
            "jaccard": 1,
            "mod_score": 1.23,
        },
        {
            "query": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "target": "9stFB1fGjCdZVZWHLVI3OD4A_DV3WcV6",
            "query_n_domains": 31,
            "target_n_domains": 1,
            "levenshtein": 0.0,
            "jaccard": 0.0,
            "mod_score": 0,
        },
        {
            "query": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "target": "MSHRSCZfdBJP8vdJdaXfeZrThH_4EUMm",
            "query_n_domains": 31,
            "target_n_domains": 2,
            "levenshtein": 0.0,
            "jaccard": 0.0,
            "mod_score": 0,
        },
        {
            "query": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "target": "Qi23auOUTcBzWTmDHuGinrzIuqH7-zVn",
            "query_n_domains": 31,
            "target_n_domains": 1,
            "levenshtein": 0.0,
            "jaccard": 0.0,
            "mod_score": 0,
        },
        {
            "query": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "target": "nk3UJUyLaWr8LochtohJ9L6eugdChZL9",
            "query_n_domains": 31,
            "target_n_domains": 0,
            "levenshtein": 0,
            "jaccard": 0,
            "mod_score": 0,
        },
        {
            "query": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "target": "du1Ncfm5UYiFYgDWD8KW1AQJNHlAcVXL",
            "query_n_domains": 31,
            "target_n_domains": 1,
            "levenshtein": 0.0,
            "jaccard": 0.0,
            "mod_score": 0,
        },
        {
            "query": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "target": "5IYMhENey2WCMrKPUz3AqBIZuFSv6DPP",
            "query_n_domains": 31,
            "target_n_domains": 0,
            "levenshtein": 0,
            "jaccard": 0,
            "mod_score": 0,
        },
        {
            "query": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "target": "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
            "query_n_domains": 31,
            "target_n_domains": 2,
            "levenshtein": 0.06,
            "jaccard": 0.2,
            "mod_score": 0.16,
        },
        {
            "query": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "target": "WbViYzQw8y-XfCQMgQXkedGduNMJPa14",
            "query_n_domains": 31,
            "target_n_domains": 3,
            "levenshtein": 0.0,
            "jaccard": 0.0,
            "mod_score": 0,
        },
        {
            "query": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "target": "4_8182J88axMDpFJBZI6kLNJAu8Ittm3",
            "query_n_domains": 31,
            "target_n_domains": 2,
            "levenshtein": 0.0,
            "jaccard": 0.0,
            "mod_score": 0,
        },
        {
            "query": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "target": "IRqRpDzrGB9UhHJD6AzDq_6Xupj00Nte",
            "query_n_domains": 31,
            "target_n_domains": 1,
            "levenshtein": 0.0,
            "jaccard": 0.0,
            "mod_score": 0,
        },
        {
            "query": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "target": "ewjVum5PbpEJA4rl-BfnCAypKl5HXb7x",
            "query_n_domains": 31,
            "target_n_domains": 0,
            "levenshtein": 0,
            "jaccard": 0,
            "mod_score": 0,
        },
    ]
    assert [i.to_dict() for i in gen] == expected


def test_CompareDomains_compare_many_to_many():
    sg_object = SocialGene()
    sg_object.parse(gbk_path)
    protein_id_list = list(sg_object.proteins.keys())
    sg_object.annotate_proteins_with_hmmscan(
        protein_id_list=protein_id_list, hmm_filepath=hmm_path, cpus=1
    )
    a = CompareDomains()
    p1 = list(sg_object.proteins.values())
    p2 = list(sg_object.proteins.values())
    gen = a.compare_many_to_many(p1, p2)
    expected = [
        {
            "query": "Tdc2m3PRLsyEzjwyux6BF4arDy2mQ_Bl",
            "target": "Tdc2m3PRLsyEzjwyux6BF4arDy2mQ_Bl",
            "query_n_domains": 2,
            "target_n_domains": 2,
            "levenshtein": 1,
            "jaccard": 1,
            "mod_score": 1.5,
        },
        {
            "query": "ptq1NGhBcUp3TIEqvAUxnnp4LOKwINvn",
            "target": "ptq1NGhBcUp3TIEqvAUxnnp4LOKwINvn",
            "query_n_domains": 0,
            "target_n_domains": 0,
            "levenshtein": 1,
            "jaccard": 1,
            "mod_score": 1.5,
        },
        {
            "query": "-l7xLyFZbiZENPLq_GML8JyTRF1Srawr",
            "target": "-l7xLyFZbiZENPLq_GML8JyTRF1Srawr",
            "query_n_domains": 0,
            "target_n_domains": 0,
            "levenshtein": 1,
            "jaccard": 1,
            "mod_score": 1.5,
        },
        {
            "query": "T_DzOorDp3ROhRRBtuXP3xyAPorpTVD0",
            "target": "T_DzOorDp3ROhRRBtuXP3xyAPorpTVD0",
            "query_n_domains": 0,
            "target_n_domains": 0,
            "levenshtein": 1,
            "jaccard": 1,
            "mod_score": 1.5,
        },
        {
            "query": "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
            "target": "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
            "query_n_domains": 14,
            "target_n_domains": 14,
            "levenshtein": 1,
            "jaccard": 1,
            "mod_score": 1.5,
        },
        {
            "query": "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
            "target": "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
            "query_n_domains": 14,
            "target_n_domains": 27,
            "levenshtein": 0.41,
            "jaccard": 0.73,
            "mod_score": 0.77,
        },
        {
            "query": "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
            "target": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "query_n_domains": 14,
            "target_n_domains": 31,
            "levenshtein": 0.35,
            "jaccard": 0.46,
            "mod_score": 0.59,
        },
        {
            "query": "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
            "target": "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
            "query_n_domains": 14,
            "target_n_domains": 7,
            "levenshtein": 0.36,
            "jaccard": 0.5,
            "mod_score": 0.61,
        },
        {
            "query": "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
            "target": "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
            "query_n_domains": 14,
            "target_n_domains": 40,
            "levenshtein": 0.28,
            "jaccard": 0.46,
            "mod_score": 0.51,
        },
        {
            "query": "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
            "target": "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
            "query_n_domains": 14,
            "target_n_domains": 2,
            "levenshtein": 0.07,
            "jaccard": 0.1,
            "mod_score": 0.12,
        },
        {
            "query": "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
            "target": "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
            "query_n_domains": 27,
            "target_n_domains": 14,
            "levenshtein": 0.41,
            "jaccard": 0.73,
            "mod_score": 0.77,
        },
        {
            "query": "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
            "target": "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
            "query_n_domains": 27,
            "target_n_domains": 27,
            "levenshtein": 1,
            "jaccard": 1,
            "mod_score": 1.5,
        },
        {
            "query": "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
            "target": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "query_n_domains": 27,
            "target_n_domains": 31,
            "levenshtein": 0.68,
            "jaccard": 0.67,
            "mod_score": 1.01,
        },
        {
            "query": "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
            "target": "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
            "query_n_domains": 27,
            "target_n_domains": 7,
            "levenshtein": 0.26,
            "jaccard": 0.6,
            "mod_score": 0.56,
        },
        {
            "query": "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
            "target": "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
            "query_n_domains": 27,
            "target_n_domains": 40,
            "levenshtein": 0.57,
            "jaccard": 0.67,
            "mod_score": 0.91,
        },
        {
            "query": "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
            "target": "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
            "query_n_domains": 27,
            "target_n_domains": 2,
            "levenshtein": 0.07,
            "jaccard": 0.2,
            "mod_score": 0.17,
        },
        {
            "query": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "target": "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
            "query_n_domains": 31,
            "target_n_domains": 14,
            "levenshtein": 0.35,
            "jaccard": 0.46,
            "mod_score": 0.59,
        },
        {
            "query": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "target": "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
            "query_n_domains": 31,
            "target_n_domains": 27,
            "levenshtein": 0.68,
            "jaccard": 0.67,
            "mod_score": 1.01,
        },
        {
            "query": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "target": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "query_n_domains": 31,
            "target_n_domains": 31,
            "levenshtein": 1,
            "jaccard": 1,
            "mod_score": 1.5,
        },
        {
            "query": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "target": "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
            "query_n_domains": 31,
            "target_n_domains": 7,
            "levenshtein": 0.23,
            "jaccard": 0.6,
            "mod_score": 0.53,
        },
        {
            "query": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "target": "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
            "query_n_domains": 31,
            "target_n_domains": 40,
            "levenshtein": 0.72,
            "jaccard": 1,
            "mod_score": 1.23,
        },
        {
            "query": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "target": "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
            "query_n_domains": 31,
            "target_n_domains": 2,
            "levenshtein": 0.06,
            "jaccard": 0.2,
            "mod_score": 0.16,
        },
        {
            "query": "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
            "target": "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
            "query_n_domains": 7,
            "target_n_domains": 14,
            "levenshtein": 0.36,
            "jaccard": 0.5,
            "mod_score": 0.61,
        },
        {
            "query": "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
            "target": "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
            "query_n_domains": 7,
            "target_n_domains": 27,
            "levenshtein": 0.26,
            "jaccard": 0.6,
            "mod_score": 0.56,
        },
        {
            "query": "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
            "target": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "query_n_domains": 7,
            "target_n_domains": 31,
            "levenshtein": 0.23,
            "jaccard": 0.6,
            "mod_score": 0.53,
        },
        {
            "query": "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
            "target": "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
            "query_n_domains": 7,
            "target_n_domains": 7,
            "levenshtein": 1,
            "jaccard": 1,
            "mod_score": 1.5,
        },
        {
            "query": "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
            "target": "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
            "query_n_domains": 7,
            "target_n_domains": 40,
            "levenshtein": 0.18,
            "jaccard": 0.6,
            "mod_score": 0.48,
        },
        {
            "query": "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
            "target": "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
            "query_n_domains": 7,
            "target_n_domains": 2,
            "levenshtein": 0.14,
            "jaccard": 0.33,
            "mod_score": 0.31,
        },
        {
            "query": "DTee9G4M8sEfnM4HaPfI37rT74pq7M_G",
            "target": "DTee9G4M8sEfnM4HaPfI37rT74pq7M_G",
            "query_n_domains": 1,
            "target_n_domains": 1,
            "levenshtein": 1,
            "jaccard": 1,
            "mod_score": 1.5,
        },
        {
            "query": "DTee9G4M8sEfnM4HaPfI37rT74pq7M_G",
            "target": "Qi23auOUTcBzWTmDHuGinrzIuqH7-zVn",
            "query_n_domains": 1,
            "target_n_domains": 1,
            "levenshtein": 1.0,
            "jaccard": 1,
            "mod_score": 1.5,
        },
        {
            "query": "mB22-i4RqtslyO7_HappM4rJ4Z2Qbkfn",
            "target": "mB22-i4RqtslyO7_HappM4rJ4Z2Qbkfn",
            "query_n_domains": 1,
            "target_n_domains": 1,
            "levenshtein": 1,
            "jaccard": 1,
            "mod_score": 1.5,
        },
        {
            "query": "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
            "target": "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
            "query_n_domains": 40,
            "target_n_domains": 14,
            "levenshtein": 0.28,
            "jaccard": 0.46,
            "mod_score": 0.51,
        },
        {
            "query": "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
            "target": "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
            "query_n_domains": 40,
            "target_n_domains": 27,
            "levenshtein": 0.57,
            "jaccard": 0.67,
            "mod_score": 0.91,
        },
        {
            "query": "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
            "target": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "query_n_domains": 40,
            "target_n_domains": 31,
            "levenshtein": 0.72,
            "jaccard": 1,
            "mod_score": 1.23,
        },
        {
            "query": "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
            "target": "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
            "query_n_domains": 40,
            "target_n_domains": 7,
            "levenshtein": 0.18,
            "jaccard": 0.6,
            "mod_score": 0.48,
        },
        {
            "query": "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
            "target": "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
            "query_n_domains": 40,
            "target_n_domains": 40,
            "levenshtein": 1,
            "jaccard": 1,
            "mod_score": 1.5,
        },
        {
            "query": "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
            "target": "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
            "query_n_domains": 40,
            "target_n_domains": 2,
            "levenshtein": 0.05,
            "jaccard": 0.2,
            "mod_score": 0.15,
        },
        {
            "query": "9stFB1fGjCdZVZWHLVI3OD4A_DV3WcV6",
            "target": "9stFB1fGjCdZVZWHLVI3OD4A_DV3WcV6",
            "query_n_domains": 1,
            "target_n_domains": 1,
            "levenshtein": 1,
            "jaccard": 1,
            "mod_score": 1.5,
        },
        {
            "query": "MSHRSCZfdBJP8vdJdaXfeZrThH_4EUMm",
            "target": "MSHRSCZfdBJP8vdJdaXfeZrThH_4EUMm",
            "query_n_domains": 2,
            "target_n_domains": 2,
            "levenshtein": 1,
            "jaccard": 1,
            "mod_score": 1.5,
        },
        {
            "query": "Qi23auOUTcBzWTmDHuGinrzIuqH7-zVn",
            "target": "DTee9G4M8sEfnM4HaPfI37rT74pq7M_G",
            "query_n_domains": 1,
            "target_n_domains": 1,
            "levenshtein": 1.0,
            "jaccard": 1,
            "mod_score": 1.5,
        },
        {
            "query": "Qi23auOUTcBzWTmDHuGinrzIuqH7-zVn",
            "target": "Qi23auOUTcBzWTmDHuGinrzIuqH7-zVn",
            "query_n_domains": 1,
            "target_n_domains": 1,
            "levenshtein": 1,
            "jaccard": 1,
            "mod_score": 1.5,
        },
        {
            "query": "nk3UJUyLaWr8LochtohJ9L6eugdChZL9",
            "target": "nk3UJUyLaWr8LochtohJ9L6eugdChZL9",
            "query_n_domains": 0,
            "target_n_domains": 0,
            "levenshtein": 1,
            "jaccard": 1,
            "mod_score": 1.5,
        },
        {
            "query": "du1Ncfm5UYiFYgDWD8KW1AQJNHlAcVXL",
            "target": "du1Ncfm5UYiFYgDWD8KW1AQJNHlAcVXL",
            "query_n_domains": 1,
            "target_n_domains": 1,
            "levenshtein": 1,
            "jaccard": 1,
            "mod_score": 1.5,
        },
        {
            "query": "5IYMhENey2WCMrKPUz3AqBIZuFSv6DPP",
            "target": "5IYMhENey2WCMrKPUz3AqBIZuFSv6DPP",
            "query_n_domains": 0,
            "target_n_domains": 0,
            "levenshtein": 1,
            "jaccard": 1,
            "mod_score": 1.5,
        },
        {
            "query": "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
            "target": "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
            "query_n_domains": 2,
            "target_n_domains": 14,
            "levenshtein": 0.07,
            "jaccard": 0.1,
            "mod_score": 0.12,
        },
        {
            "query": "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
            "target": "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
            "query_n_domains": 2,
            "target_n_domains": 27,
            "levenshtein": 0.07,
            "jaccard": 0.2,
            "mod_score": 0.17,
        },
        {
            "query": "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
            "target": "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "query_n_domains": 2,
            "target_n_domains": 31,
            "levenshtein": 0.06,
            "jaccard": 0.2,
            "mod_score": 0.16,
        },
        {
            "query": "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
            "target": "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
            "query_n_domains": 2,
            "target_n_domains": 7,
            "levenshtein": 0.14,
            "jaccard": 0.33,
            "mod_score": 0.31,
        },
        {
            "query": "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
            "target": "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
            "query_n_domains": 2,
            "target_n_domains": 40,
            "levenshtein": 0.05,
            "jaccard": 0.2,
            "mod_score": 0.15,
        },
        {
            "query": "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
            "target": "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
            "query_n_domains": 2,
            "target_n_domains": 2,
            "levenshtein": 1,
            "jaccard": 1,
            "mod_score": 1.5,
        },
        {
            "query": "WbViYzQw8y-XfCQMgQXkedGduNMJPa14",
            "target": "WbViYzQw8y-XfCQMgQXkedGduNMJPa14",
            "query_n_domains": 3,
            "target_n_domains": 3,
            "levenshtein": 1,
            "jaccard": 1,
            "mod_score": 1.5,
        },
        {
            "query": "4_8182J88axMDpFJBZI6kLNJAu8Ittm3",
            "target": "4_8182J88axMDpFJBZI6kLNJAu8Ittm3",
            "query_n_domains": 2,
            "target_n_domains": 2,
            "levenshtein": 1,
            "jaccard": 1,
            "mod_score": 1.5,
        },
        {
            "query": "IRqRpDzrGB9UhHJD6AzDq_6Xupj00Nte",
            "target": "IRqRpDzrGB9UhHJD6AzDq_6Xupj00Nte",
            "query_n_domains": 1,
            "target_n_domains": 1,
            "levenshtein": 1,
            "jaccard": 1,
            "mod_score": 1.5,
        },
        {
            "query": "ewjVum5PbpEJA4rl-BfnCAypKl5HXb7x",
            "target": "ewjVum5PbpEJA4rl-BfnCAypKl5HXb7x",
            "query_n_domains": 0,
            "target_n_domains": 0,
            "levenshtein": 1,
            "jaccard": 1,
            "mod_score": 1.5,
        },
    ]

    assert [i.to_dict() for i in gen] == expected
    gen = a.compare_many_to_many(p1, p2, filter_non_hits=False)
    # 22 inputs * 22 inputs = 484
    assert len(list(gen)) == 484


def test_CompareDomains_compare_all_to_all_parallel():
    sg_object = SocialGene()
    sg_object.parse(gbk_path)
    protein_id_list = list(sg_object.proteins.keys())
    sg_object.annotate_proteins_with_hmmscan(
        protein_id_list=protein_id_list, hmm_filepath=hmm_path, cpus=1
    )
    a = CompareDomains()
    p2 = list(sg_object.proteins.values())
    a.compare_all_to_all_parallel(p2)
