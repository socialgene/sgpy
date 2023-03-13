import os
from socialgene.base.socialgene import SocialGene

FIXTURE_DIR = os.path.dirname(os.path.realpath(__file__))
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)
FIXTURE_DIR = os.path.join(FIXTURE_DIR, "data")
FIXTURE_DIR = os.path.join(FIXTURE_DIR, "lagriamide_mibig_bgc0001946.gbk")


def test_protein_info_table_1():
    expected_protein_info_table = [
        ("Tdc2m3PRLsyEzjwyux6BF4arDy2mQ_Bl", "AXA20086.1", "sigma-70 RpoE", 174),
        (
            "ptq1NGhBcUp3TIEqvAUxnnp4LOKwINvn",
            "AXA20087.1",
            "competence protein ComEC",
            64,
        ),
        ("-l7xLyFZbiZENPLq_GML8JyTRF1Srawr", "AXA20088.1", "transposase", 36),
        ("T_DzOorDp3ROhRRBtuXP3xyAPorpTVD0", "AXA20089.1", "hypothetical protein", 96),
        (
            "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
            "AXA20090.1",
            "hybrid trans-AT PKS/NRPS LgaA",
            3553,
        ),
        (
            "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
            "AXA20091.1",
            "hybrid trans-AT PKS/NRPS LgaB",
            6405,
        ),
        ("Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5", "AXA20092.1", "trans-AT PKS LgaC", 6799),
        ("RyDIaUZc_b21_kQalx7J3yNO4l5f-439", "AXA20093.1", "trans-AT PKS LgaD", 1279),
        ("DTee9G4M8sEfnM4HaPfI37rT74pq7M_G", "AXA20094.1", "acyltransferase LgaE", 370),
        ("mB22-i4RqtslyO7_HappM4rJ4Z2Qbkfn", "AXA20095.1", "enoylreductase LgaF", 461),
        ("iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh", "AXA20096.1", "trans-AT PKS LgaG", 8904),
        (
            "9stFB1fGjCdZVZWHLVI3OD4A_DV3WcV6",
            "AXA20097.1",
            "nuclear transport factor 2 (NTF2)-like protein LgaL",
            131,
        ),
        (
            "MSHRSCZfdBJP8vdJdaXfeZrThH_4EUMm",
            "AXA20098.1",
            "MATE family efflux transporter LgaH",
            457,
        ),
        ("Qi23auOUTcBzWTmDHuGinrzIuqH7-zVn", "AXA20099.1", "acylhydrolase LgaI", 337),
        ("nk3UJUyLaWr8LochtohJ9L6eugdChZL9", "AXA20100.1", "transposase", 128),
        ("du1Ncfm5UYiFYgDWD8KW1AQJNHlAcVXL", "AXA20101.1", "cytochrome P450 LgaJ", 394),
        ("5IYMhENey2WCMrKPUz3AqBIZuFSv6DPP", "AXA20102.1", "LgaM", 244),
        ("iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L", "AXA20103.1", "ketoreductase LgaK", 273),
        (
            "WbViYzQw8y-XfCQMgQXkedGduNMJPa14",
            "AXA20104.1",
            "phosphoenolpyruvate-protein phosphotransferase",
            587,
        ),
        (
            "4_8182J88axMDpFJBZI6kLNJAu8Ittm3",
            "AXA20105.1",
            "squalene--hopene cyclase",
            679,
        ),
        (
            "IRqRpDzrGB9UhHJD6AzDq_6Xupj00Nte",
            "AXA20106.1",
            "all-trans-phytoene synthase",
            281,
        ),
        ("ewjVum5PbpEJA4rl-BfnCAypKl5HXb7x", "AXA20107.1", "hypothetical protein", 155),
    ]
    sg_obj = SocialGene()
    sg_obj.parse(FIXTURE_DIR)
    sg_obj.proteins
    assert [i for i in sg_obj.protein_info_table()] == expected_protein_info_table


def test_protein_info_table_2():
    expected_protein_info_table = [
        ("Tdc2m3PRLsyEzjwyux6BF4arDy2mQ_Bl", "AXA20086.1", "sigma-70 RpoE", None),
        (
            "ptq1NGhBcUp3TIEqvAUxnnp4LOKwINvn",
            "AXA20087.1",
            "competence protein ComEC",
            64,
        ),
        ("-l7xLyFZbiZENPLq_GML8JyTRF1Srawr", "AXA20088.1", "transposase", 36),
        ("T_DzOorDp3ROhRRBtuXP3xyAPorpTVD0", "AXA20089.1", "hypothetical protein", 96),
        (
            "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
            "AXA20090.1",
            "hybrid trans-AT PKS/NRPS LgaA",
            3553,
        ),
        (
            "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
            "AXA20091.1",
            "hybrid trans-AT PKS/NRPS LgaB",
            6405,
        ),
        ("Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5", "AXA20092.1", "trans-AT PKS LgaC", 6799),
        ("RyDIaUZc_b21_kQalx7J3yNO4l5f-439", "AXA20093.1", "trans-AT PKS LgaD", None),
        ("DTee9G4M8sEfnM4HaPfI37rT74pq7M_G", "AXA20094.1", "acyltransferase LgaE", 370),
        ("mB22-i4RqtslyO7_HappM4rJ4Z2Qbkfn", "AXA20095.1", "enoylreductase LgaF", 461),
        ("iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh", "AXA20096.1", "trans-AT PKS LgaG", 8904),
        (
            "9stFB1fGjCdZVZWHLVI3OD4A_DV3WcV6",
            "AXA20097.1",
            "nuclear transport factor 2 (NTF2)-like protein LgaL",
            131,
        ),
        (
            "MSHRSCZfdBJP8vdJdaXfeZrThH_4EUMm",
            "AXA20098.1",
            "MATE family efflux transporter LgaH",
            457,
        ),
        ("Qi23auOUTcBzWTmDHuGinrzIuqH7-zVn", "AXA20099.1", "acylhydrolase LgaI", 337),
        ("nk3UJUyLaWr8LochtohJ9L6eugdChZL9", "AXA20100.1", "transposase", 128),
        ("du1Ncfm5UYiFYgDWD8KW1AQJNHlAcVXL", "AXA20101.1", "cytochrome P450 LgaJ", 394),
        ("5IYMhENey2WCMrKPUz3AqBIZuFSv6DPP", "AXA20102.1", "LgaM", 244),
        ("iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L", "AXA20103.1", "ketoreductase LgaK", 273),
        (
            "WbViYzQw8y-XfCQMgQXkedGduNMJPa14",
            "AXA20104.1",
            "phosphoenolpyruvate-protein phosphotransferase",
            587,
        ),
        (
            "4_8182J88axMDpFJBZI6kLNJAu8Ittm3",
            "AXA20105.1",
            "squalene--hopene cyclase",
            679,
        ),
        (
            "IRqRpDzrGB9UhHJD6AzDq_6Xupj00Nte",
            "AXA20106.1",
            "all-trans-phytoene synthase",
            281,
        ),
        ("ewjVum5PbpEJA4rl-BfnCAypKl5HXb7x", "AXA20107.1", "hypothetical protein", 155),
    ]
    sg_obj = SocialGene()
    sg_obj.parse(FIXTURE_DIR)
    sg_obj.proteins["Tdc2m3PRLsyEzjwyux6BF4arDy2mQ_Bl"].sequence = None
    sg_obj.proteins["RyDIaUZc_b21_kQalx7J3yNO4l5f-439"].sequence = None
    assert [i for i in sg_obj.protein_info_table()] == expected_protein_info_table
