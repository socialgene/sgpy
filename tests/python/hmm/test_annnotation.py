import os

from socialgene.base.socialgene import SocialGene
from socialgene.hmm.hmmer import HMMER

FIXTURE_DIR = os.path.dirname(os.path.realpath(__file__))
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)
FIXTURE_DIR = os.path.join(FIXTURE_DIR, "data", "test_genomes")
gbk_path = os.path.join(FIXTURE_DIR, "lagriamide_mibig_bgc0001946.gbk")


FIXTURE_DIR = os.path.dirname(os.path.realpath(__file__))
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)
FIXTURE_DIR = os.path.join(FIXTURE_DIR, "data", "hmms")
hmm_path = os.path.join(FIXTURE_DIR, "pks.hmm")


expected_proteins = [
    "-l7xLyFZbiZENPLq_GML8JyTRF1Srawr",
    "4_8182J88axMDpFJBZI6kLNJAu8Ittm3",
    "5IYMhENey2WCMrKPUz3AqBIZuFSv6DPP",
    "9stFB1fGjCdZVZWHLVI3OD4A_DV3WcV6",
    "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
    "DTee9G4M8sEfnM4HaPfI37rT74pq7M_G",
    "IRqRpDzrGB9UhHJD6AzDq_6Xupj00Nte",
    "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
    "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
    "MSHRSCZfdBJP8vdJdaXfeZrThH_4EUMm",
    "Qi23auOUTcBzWTmDHuGinrzIuqH7-zVn",
    "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
    "T_DzOorDp3ROhRRBtuXP3xyAPorpTVD0",
    "Tdc2m3PRLsyEzjwyux6BF4arDy2mQ_Bl",
    "WbViYzQw8y-XfCQMgQXkedGduNMJPa14",
    "du1Ncfm5UYiFYgDWD8KW1AQJNHlAcVXL",
    "ewjVum5PbpEJA4rl-BfnCAypKl5HXb7x",
    "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
    "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
    "mB22-i4RqtslyO7_HappM4rJ4Z2Qbkfn",
    "nk3UJUyLaWr8LochtohJ9L6eugdChZL9",
    "ptq1NGhBcUp3TIEqvAUxnnp4LOKwINvn",
]


def test_hmmscan():
    sg_object = SocialGene()
    sg_object.parse(gbk_path)
    protein_id_list = list(sg_object.proteins.keys())

    h = HMMER()
    h.hmmpress(hmm_path)
    sg_object.annotate_proteins_with_hmmscan(
        protein_id_list=protein_id_list, hmm_directory=os.path.dirname(hmm_path), cpus=1
    )
    prots = list(sg_object.proteins.keys())
    prots.sort()
    assert prots == expected_proteins
    assert sorted(
        sg_object.proteins["Tdc2m3PRLsyEzjwyux6BF4arDy2mQ_Bl"].domain_vector
    ) == sorted(
        [
            "DIJQMpAiLKGDPgcpc1IuBzFdf7FhTYu5",
            "DIJQMpAiLKGDPgcpc1IuBzFdf7FhTYu5",
            "xJwofaGb0EIZrSxSeZL5xS6thEM7ck7U",
            "xJwofaGb0EIZrSxSeZL5xS6thEM7ck7U",
        ]
    )


def test_hmmscan2():
    sg_object = SocialGene()
    sg_object.parse(gbk_path)
    protein_id_list = list(sg_object.proteins.keys())
    h = HMMER()
    h.hmmpress(hmm_path)
    sg_object.annotate_proteins_with_hmmscan(
        protein_id_list=protein_id_list, hmm_directory=os.path.dirname(hmm_path), cpus=1
    )
    assert [
        dict(i.__dict__)
        for i in sg_object.proteins[
            "Tdc2m3PRLsyEzjwyux6BF4arDy2mQ_Bl"
        ].domain_list_sorted_by_mean_envelope_position
    ] == [
        {
            "hmm_id": "xJwofaGb0EIZrSxSeZL5xS6thEM7ck7U",
            "env_from": 2,
            "env_to": 11,
            "seq_pro_score": 45.5,
            "evalue": -8,
            "i_evalue": 7,
            "domain_bias": 0.3,
            "domain_score": -3.3,
            "seq_pro_bias": 0.1,
            "hmm_from": 46,
            "hmm_to": 54,
            "ali_from": 3,
            "ali_to": 11,
            "exponentialized": True,
        },
        {
            "hmm_id": "DIJQMpAiLKGDPgcpc1IuBzFdf7FhTYu5",
            "env_from": 4,
            "env_to": 66,
            "seq_pro_score": 28.6,
            "evalue": -3,
            "i_evalue": -3,
            "domain_bias": 0.1,
            "domain_score": 26.9,
            "seq_pro_bias": 0.2,
            "hmm_from": 19,
            "hmm_to": 69,
            "ali_from": 16,
            "ali_to": 65,
            "exponentialized": True,
        },
        {
            "hmm_id": "DIJQMpAiLKGDPgcpc1IuBzFdf7FhTYu5",
            "env_from": 93,
            "env_to": 128,
            "seq_pro_score": 28.6,
            "evalue": -3,
            "i_evalue": 7,
            "domain_bias": 0.0,
            "domain_score": -2.2,
            "seq_pro_bias": 0.2,
            "hmm_from": 38,
            "hmm_to": 62,
            "ali_from": 98,
            "ali_to": 122,
            "exponentialized": True,
        },
        {
            "hmm_id": "xJwofaGb0EIZrSxSeZL5xS6thEM7ck7U",
            "env_from": 95,
            "env_to": 147,
            "seq_pro_score": 45.5,
            "evalue": -8,
            "i_evalue": -8,
            "domain_bias": 0.1,
            "domain_score": 45.5,
            "seq_pro_bias": 0.1,
            "hmm_from": 3,
            "hmm_to": 54,
            "ali_from": 96,
            "ali_to": 147,
            "exponentialized": True,
        },
    ]
