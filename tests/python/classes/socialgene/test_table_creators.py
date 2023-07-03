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
            "37C2BC6B34544B2E",
            "AXA20086.1",
            "sigma-70 RpoE",
            174,
        ),
        (
            "ptq1NGhBcUp3TIEqvAUxnnp4LOKwINvn",
            "24885C21E4F10395",
            "AXA20087.1",
            "competence protein ComEC",
            64,
        ),
        (
            "-l7xLyFZbiZENPLq_GML8JyTRF1Srawr",
            "5DB7DF559DF798B9",
            "AXA20088.1",
            "transposase",
            36,
        ),
        (
            "T_DzOorDp3ROhRRBtuXP3xyAPorpTVD0",
            "380F00AA06102D54",
            "AXA20089.1",
            "hypothetical protein",
            96,
        ),
        (
            "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
            "1CF65176FCF98848",
            "AXA20090.1",
            "hybrid trans-AT PKS/NRPS LgaA",
            3553,
        ),
        (
            "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
            "3BFF8AA6C1E11053",
            "AXA20091.1",
            "hybrid trans-AT PKS/NRPS LgaB",
            6405,
        ),
        (
            "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "5A6503F0ADBBF5DD",
            "AXA20092.1",
            "trans-AT PKS LgaC",
            6799,
        ),
        (
            "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
            "B5363BCBA49EA647",
            "AXA20093.1",
            "trans-AT PKS LgaD",
            1279,
        ),
        (
            "DTee9G4M8sEfnM4HaPfI37rT74pq7M_G",
            "93E706B4F7EB53AF",
            "AXA20094.1",
            "acyltransferase LgaE",
            370,
        ),
        (
            "mB22-i4RqtslyO7_HappM4rJ4Z2Qbkfn",
            "0FB2B60341C61380",
            "AXA20095.1",
            "enoylreductase LgaF",
            461,
        ),
        (
            "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
            "53C03058C526F0D2",
            "AXA20096.1",
            "trans-AT PKS LgaG",
            8904,
        ),
        (
            "9stFB1fGjCdZVZWHLVI3OD4A_DV3WcV6",
            "64F56D5193D6A584",
            "AXA20097.1",
            "nuclear transport factor 2 (NTF2)-like protein LgaL",
            131,
        ),
        (
            "MSHRSCZfdBJP8vdJdaXfeZrThH_4EUMm",
            "A4115D03F2A9920B",
            "AXA20098.1",
            "MATE family efflux transporter LgaH",
            457,
        ),
        (
            "Qi23auOUTcBzWTmDHuGinrzIuqH7-zVn",
            "631C4C8531E0AE0F",
            "AXA20099.1",
            "acylhydrolase LgaI",
            337,
        ),
        (
            "nk3UJUyLaWr8LochtohJ9L6eugdChZL9",
            "00CD39A87B4E5579",
            "AXA20100.1",
            "transposase",
            128,
        ),
        (
            "du1Ncfm5UYiFYgDWD8KW1AQJNHlAcVXL",
            "D9B06B39B6ED339A",
            "AXA20101.1",
            "cytochrome P450 LgaJ",
            394,
        ),
        (
            "5IYMhENey2WCMrKPUz3AqBIZuFSv6DPP",
            "FF58B5837B0F03BE",
            "AXA20102.1",
            "LgaM",
            244,
        ),
        (
            "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
            "C16214B7E088BFA5",
            "AXA20103.1",
            "ketoreductase LgaK",
            273,
        ),
        (
            "WbViYzQw8y-XfCQMgQXkedGduNMJPa14",
            "2CC09AA68B15064E",
            "AXA20104.1",
            "phosphoenolpyruvate-protein phosphotransferase",
            587,
        ),
        (
            "4_8182J88axMDpFJBZI6kLNJAu8Ittm3",
            "F0D86567887A4777",
            "AXA20105.1",
            "squalene--hopene cyclase",
            679,
        ),
        (
            "IRqRpDzrGB9UhHJD6AzDq_6Xupj00Nte",
            "5CD80CF74497D50C",
            "AXA20106.1",
            "all-trans-phytoene synthase",
            281,
        ),
        (
            "ewjVum5PbpEJA4rl-BfnCAypKl5HXb7x",
            "E3DBD3CB4B0B011F",
            "AXA20107.1",
            "hypothetical protein",
            155,
        ),
    ]
    sg_obj = SocialGene()
    sg_obj.parse(FIXTURE_DIR)
    sg_obj.proteins
    assert [i for i in sg_obj.protein_info_table()] == expected_protein_info_table
