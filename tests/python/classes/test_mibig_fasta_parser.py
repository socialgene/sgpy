import os
import tempfile

from socialgene.base.socialgene import SocialGene

from .test_mibig_gbk_parser import PROTEIN_DICT

DIRECTORY_OF_THIS_FILE = os.path.dirname(os.path.realpath(__file__))
FIXTURE_DIR = os.path.dirname(DIRECTORY_OF_THIS_FILE)

FIXTURE_DIR = os.path.join(FIXTURE_DIR, "data", "test_genomes")


def test_fasta_file_parse():
    # just create the fasta from the genbank file
    with tempfile.NamedTemporaryFile(suffix=".faa") as fp:
        sg_object = SocialGene()
        gbk_path = os.path.join(FIXTURE_DIR, "lagriamide_mibig_bgc0001946.gbk")
        sg_object.parse(gbk_path)
        sg_object.write_fasta(outpath=fp.name)
        fasta_object = SocialGene()
        fasta_object.parse(fp.name)

    protein_parse_results = {
        k: [v.description, v.external_protein_id, v.domains]
        for k, v in fasta_object.proteins.items()
    }
    assert protein_parse_results == {
        "Tdc2m3PRLsyEzjwyux6BF4arDy2mQ_Bl": [
            "Tdc2m3PRLsyEzjwyux6BF4arDy2mQ_Bl",
            "Tdc2m3PRLsyEzjwyux6BF4arDy2mQ_Bl",
            set(),
        ],
        "ptq1NGhBcUp3TIEqvAUxnnp4LOKwINvn": [
            "ptq1NGhBcUp3TIEqvAUxnnp4LOKwINvn",
            "ptq1NGhBcUp3TIEqvAUxnnp4LOKwINvn",
            set(),
        ],
        "-l7xLyFZbiZENPLq_GML8JyTRF1Srawr": [
            "-l7xLyFZbiZENPLq_GML8JyTRF1Srawr",
            "-l7xLyFZbiZENPLq_GML8JyTRF1Srawr",
            set(),
        ],
        "T_DzOorDp3ROhRRBtuXP3xyAPorpTVD0": [
            "T_DzOorDp3ROhRRBtuXP3xyAPorpTVD0",
            "T_DzOorDp3ROhRRBtuXP3xyAPorpTVD0",
            set(),
        ],
        "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl": [
            "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
            "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl",
            set(),
        ],
        "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ": [
            "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
            "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ",
            set(),
        ],
        "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5": [
            "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5",
            set(),
        ],
        "RyDIaUZc_b21_kQalx7J3yNO4l5f-439": [
            "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
            "RyDIaUZc_b21_kQalx7J3yNO4l5f-439",
            set(),
        ],
        "DTee9G4M8sEfnM4HaPfI37rT74pq7M_G": [
            "DTee9G4M8sEfnM4HaPfI37rT74pq7M_G",
            "DTee9G4M8sEfnM4HaPfI37rT74pq7M_G",
            set(),
        ],
        "mB22-i4RqtslyO7_HappM4rJ4Z2Qbkfn": [
            "mB22-i4RqtslyO7_HappM4rJ4Z2Qbkfn",
            "mB22-i4RqtslyO7_HappM4rJ4Z2Qbkfn",
            set(),
        ],
        "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh": [
            "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
            "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh",
            set(),
        ],
        "9stFB1fGjCdZVZWHLVI3OD4A_DV3WcV6": [
            "9stFB1fGjCdZVZWHLVI3OD4A_DV3WcV6",
            "9stFB1fGjCdZVZWHLVI3OD4A_DV3WcV6",
            set(),
        ],
        "MSHRSCZfdBJP8vdJdaXfeZrThH_4EUMm": [
            "MSHRSCZfdBJP8vdJdaXfeZrThH_4EUMm",
            "MSHRSCZfdBJP8vdJdaXfeZrThH_4EUMm",
            set(),
        ],
        "Qi23auOUTcBzWTmDHuGinrzIuqH7-zVn": [
            "Qi23auOUTcBzWTmDHuGinrzIuqH7-zVn",
            "Qi23auOUTcBzWTmDHuGinrzIuqH7-zVn",
            set(),
        ],
        "nk3UJUyLaWr8LochtohJ9L6eugdChZL9": [
            "nk3UJUyLaWr8LochtohJ9L6eugdChZL9",
            "nk3UJUyLaWr8LochtohJ9L6eugdChZL9",
            set(),
        ],
        "du1Ncfm5UYiFYgDWD8KW1AQJNHlAcVXL": [
            "du1Ncfm5UYiFYgDWD8KW1AQJNHlAcVXL",
            "du1Ncfm5UYiFYgDWD8KW1AQJNHlAcVXL",
            set(),
        ],
        "5IYMhENey2WCMrKPUz3AqBIZuFSv6DPP": [
            "5IYMhENey2WCMrKPUz3AqBIZuFSv6DPP",
            "5IYMhENey2WCMrKPUz3AqBIZuFSv6DPP",
            set(),
        ],
        "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L": [
            "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
            "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L",
            set(),
        ],
        "WbViYzQw8y-XfCQMgQXkedGduNMJPa14": [
            "WbViYzQw8y-XfCQMgQXkedGduNMJPa14",
            "WbViYzQw8y-XfCQMgQXkedGduNMJPa14",
            set(),
        ],
        "4_8182J88axMDpFJBZI6kLNJAu8Ittm3": [
            "4_8182J88axMDpFJBZI6kLNJAu8Ittm3",
            "4_8182J88axMDpFJBZI6kLNJAu8Ittm3",
            set(),
        ],
        "IRqRpDzrGB9UhHJD6AzDq_6Xupj00Nte": [
            "IRqRpDzrGB9UhHJD6AzDq_6Xupj00Nte",
            "IRqRpDzrGB9UhHJD6AzDq_6Xupj00Nte",
            set(),
        ],
        "ewjVum5PbpEJA4rl-BfnCAypKl5HXb7x": [
            "ewjVum5PbpEJA4rl-BfnCAypKl5HXb7x",
            "ewjVum5PbpEJA4rl-BfnCAypKl5HXb7x",
            set(),
        ],
    }
    assert {k: v.sequence for k, v in sg_object.proteins.items()} == PROTEIN_DICT


def test_fasta_string_parse():
    sg_object = SocialGene()
    sg_object.parse_fasta_string(">asdads\n dasfa")
    assert list(sg_object.proteins.keys())[0] == "lZA5w9NxutVBhoqjqfsX_GX_dfugQZLQ"
    protein_parse_results = [
        {k: [v.description, v.external_protein_id, v.domains]}
        for k, v in sg_object.proteins.items()
    ]
    assert protein_parse_results == [
        {"lZA5w9NxutVBhoqjqfsX_GX_dfugQZLQ": ["asdads", "asdads", set()]}
    ]
