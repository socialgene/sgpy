import os

from socialgene.base.socialgene import SocialGene

FIXTURE_DIR = os.path.dirname(os.path.realpath(__file__))
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)
FIXTURE_DIR = os.path.join(FIXTURE_DIR, "data", "test_genomes")
FIXTURE_DIR = os.path.join(FIXTURE_DIR, "lagriamide_mibig_bgc0001946.gbk")


def test_protein_info_table_1():
    expected_protein_info_table = [
        (
            "Tdc2m3PRLsyEzjwyux6BF4arDy2mQ_Bl",
            "6a8bbacb2c6d3863fc95b69bc726c0a6",
            "AXA20086.1",
            "sigma-70 RpoE",
            174,
        ),
        (
            "ptq1NGhBcUp3TIEqvAUxnnp4LOKwINvn",
            "aea834b3e35ad5ff47838864fd189429",
            "AXA20087.1",
            "competence protein ComEC",
            64,
        ),
        (
            "-l7xLyFZbiZENPLq_GML8JyTRF1Srawr",
            "3dc35df848038994ec8c560a0466d2ec",
            "AXA20088.1",
            "transposase",
            36,
        ),
        (
            "T_DzOorDp3ROhRRBtuXP3xyAPorpTVD0",
            "1d1b84699a2c9a0794f16f4bc2d175f8",
            "AXA20089.1",
            "hypothetical protein",
            96,
        ),
        (
            "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
            "f7310dfc5a871ebb64355213a5747179",
            "AXA20090.1",
            "hybrid trans-AT PKS/NRPS LgaA",
            3553,
        ),
        (
            "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
            "d10d291e7fbce5f91cfd7237c7219752",
            "AXA20091.1",
            "hybrid trans-AT PKS/NRPS LgaB",
            6405,
        ),
        (
            "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "9514359ef2bf947c0b6d916ed50c8cd5",
            "AXA20092.1",
            "trans-AT PKS LgaC",
            6799,
        ),
        (
            "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
            "6955ad1615bef63f9d21ed431a40493f",
            "AXA20093.1",
            "trans-AT PKS LgaD",
            1279,
        ),
        (
            "DTee9G4M8sEfnM4HaPfI37rT74pq7M_G",
            "05e4301e11acfcc9066e9827d984a9fe",
            "AXA20094.1",
            "acyltransferase LgaE",
            370,
        ),
        (
            "mB22-i4RqtslyO7_HappM4rJ4Z2Qbkfn",
            "13f6ea0cde2ed51885ffc359e9c4b56c",
            "AXA20095.1",
            "enoylreductase LgaF",
            461,
        ),
        (
            "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
            "9a54df1d31a5530ab510fa057c0d39a7",
            "AXA20096.1",
            "trans-AT PKS LgaG",
            8904,
        ),
        (
            "9stFB1fGjCdZVZWHLVI3OD4A_DV3WcV6",
            "6f2108956277fe58270dd62610cfcd62",
            "AXA20097.1",
            "nuclear transport factor 2 (NTF2)-like protein LgaL",
            131,
        ),
        (
            "MSHRSCZfdBJP8vdJdaXfeZrThH_4EUMm",
            "dbc45af66c60fc5e3f924333f525e728",
            "AXA20098.1",
            "MATE family efflux transporter LgaH",
            457,
        ),
        (
            "Qi23auOUTcBzWTmDHuGinrzIuqH7-zVn",
            "f9c1d6d90db6cb41128ad1f6a9b4b6cf",
            "AXA20099.1",
            "acylhydrolase LgaI",
            337,
        ),
        (
            "nk3UJUyLaWr8LochtohJ9L6eugdChZL9",
            "a694599d8141567683821798d1e4f370",
            "AXA20100.1",
            "transposase",
            128,
        ),
        (
            "du1Ncfm5UYiFYgDWD8KW1AQJNHlAcVXL",
            "454262859775ba1c026194d1b3742f2f",
            "AXA20101.1",
            "cytochrome P450 LgaJ",
            394,
        ),
        (
            "5IYMhENey2WCMrKPUz3AqBIZuFSv6DPP",
            "9ccd7dd17552cb2c1d66f4bc637eb367",
            "AXA20102.1",
            "LgaM",
            244,
        ),
        (
            "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
            "7f74bf1901415286b8e0e2dab5c31c71",
            "AXA20103.1",
            "ketoreductase LgaK",
            273,
        ),
        (
            "WbViYzQw8y-XfCQMgQXkedGduNMJPa14",
            "09d974bcebc9745c7686c350ee411fe0",
            "AXA20104.1",
            "phosphoenolpyruvate-protein phosphotransferase",
            587,
        ),
        (
            "4_8182J88axMDpFJBZI6kLNJAu8Ittm3",
            "5f888dca5594c850f88b7f7f292614d6",
            "AXA20105.1",
            "squalene--hopene cyclase",
            679,
        ),
        (
            "IRqRpDzrGB9UhHJD6AzDq_6Xupj00Nte",
            "411bd582129cbcc1434d922a176a3c31",
            "AXA20106.1",
            "all-trans-phytoene synthase",
            281,
        ),
        (
            "ewjVum5PbpEJA4rl-BfnCAypKl5HXb7x",
            "041556364317af40b5b3df5674388c53",
            "AXA20107.1",
            "hypothetical protein",
            155,
        ),
    ]
    sg_obj = SocialGene()
    sg_obj.parse(FIXTURE_DIR)
    sg_obj.proteins
    assert [i for i in sg_obj.protein_info_table()] == expected_protein_info_table
