import os
import hashlib
from socialgene.base.socialgene import SocialGene
import socialgene.base.molbio as molbio
import tempfile

FIXTURE_DIR = os.path.dirname(os.path.realpath(__file__))
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)

FIXTURE_DIR = os.path.join(FIXTURE_DIR, "data")


def test_genbank_file_correct_parse_structure():
    sg_object = SocialGene()
    gbk_path = os.path.join(FIXTURE_DIR, "lagriamide_mibig_bgc0001946.gbk")
    sg_object.parse(gbk_path)
    assert isinstance(sg_object, SocialGene)
    assert isinstance(sg_object, SocialGene)
    assert isinstance(sg_object.assemblies, dict)
    assert isinstance(sg_object.assemblies["lagriamide_mibig_bgc0001946"].loci, dict)
    # locus contains a uuid so...
    locus_key = [
        key
        for key, value in sg_object.assemblies[
            "lagriamide_mibig_bgc0001946"
        ].loci.items()
        if "bgc0001946" in key.lower()
    ][0]
    assert isinstance(
        sg_object.assemblies["lagriamide_mibig_bgc0001946"].loci[locus_key],
        molbio.Locus,
    )
    assert isinstance(sg_object.proteins, dict)
    assert isinstance(
        sg_object.proteins["Tdc2m3PRLsyEzjwyux6BF4arDy2mQ_Bl"], molbio.Protein
    )


def test_genbank_file_parse_result():
    sg_object = SocialGene()
    gbk_path = os.path.join(FIXTURE_DIR, "lagriamide_mibig_bgc0001946.gbk")
    sg_object.parse(gbk_path)
    # locus contains a uuid so...
    locus_key = [
        key
        for key, value in sg_object.assemblies[
            "lagriamide_mibig_bgc0001946"
        ].loci.items()
        if "bgc0001946" in key.lower()
    ][0]
    protein_parse_results = [
        {k: [v.description, v.other_id, v.domains]}
        for k, v in sg_object.proteins.items()
    ]
    assert protein_parse_results == [
        {
            "Tdc2m3PRLsyEzjwyux6BF4arDy2mQ_Bl": [
                "sigma-70 RpoE",
                "AXA20086.1",
                set(),
            ]
        },
        {
            "ptq1NGhBcUp3TIEqvAUxnnp4LOKwINvn": [
                "competence protein ComEC",
                "AXA20087.1",
                set(),
            ]
        },
        {
            "-l7xLyFZbiZENPLq_GML8JyTRF1Srawr": [
                "transposase",
                "AXA20088.1",
                set(),
            ]
        },
        {
            "T_DzOorDp3ROhRRBtuXP3xyAPorpTVD0": [
                "hypothetical protein",
                "AXA20089.1",
                set(),
            ]
        },
        {
            "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl": [
                "hybrid trans-AT PKS/NRPS LgaA",
                "AXA20090.1",
                set(),
            ]
        },
        {
            "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ": [
                "hybrid trans-AT PKS/NRPS LgaB",
                "AXA20091.1",
                set(),
            ]
        },
        {
            "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5": [
                "trans-AT PKS LgaC",
                "AXA20092.1",
                set(),
            ]
        },
        {
            "RyDIaUZc_b21_kQalx7J3yNO4l5f-439": [
                "trans-AT PKS LgaD",
                "AXA20093.1",
                set(),
            ]
        },
        {
            "DTee9G4M8sEfnM4HaPfI37rT74pq7M_G": [
                "acyltransferase LgaE",
                "AXA20094.1",
                set(),
            ]
        },
        {
            "mB22-i4RqtslyO7_HappM4rJ4Z2Qbkfn": [
                "enoylreductase LgaF",
                "AXA20095.1",
                set(),
            ]
        },
        {
            "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh": [
                "trans-AT PKS LgaG",
                "AXA20096.1",
                set(),
            ]
        },
        {
            "9stFB1fGjCdZVZWHLVI3OD4A_DV3WcV6": [
                "nuclear transport factor 2 (NTF2)-like protein LgaL",
                "AXA20097.1",
                set(),
            ]
        },
        {
            "MSHRSCZfdBJP8vdJdaXfeZrThH_4EUMm": [
                "MATE family efflux transporter LgaH",
                "AXA20098.1",
                set(),
            ]
        },
        {
            "Qi23auOUTcBzWTmDHuGinrzIuqH7-zVn": [
                "acylhydrolase LgaI",
                "AXA20099.1",
                set(),
            ]
        },
        {
            "nk3UJUyLaWr8LochtohJ9L6eugdChZL9": [
                "transposase",
                "AXA20100.1",
                set(),
            ]
        },
        {
            "du1Ncfm5UYiFYgDWD8KW1AQJNHlAcVXL": [
                "cytochrome P450 LgaJ",
                "AXA20101.1",
                set(),
            ]
        },
        {"5IYMhENey2WCMrKPUz3AqBIZuFSv6DPP": ["LgaM", "AXA20102.1", set()]},
        {
            "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L": [
                "ketoreductase LgaK",
                "AXA20103.1",
                set(),
            ]
        },
        {
            "WbViYzQw8y-XfCQMgQXkedGduNMJPa14": [
                "phosphoenolpyruvate-protein phosphotransferase",
                "AXA20104.1",
                set(),
            ]
        },
        {
            "4_8182J88axMDpFJBZI6kLNJAu8Ittm3": [
                "squalene--hopene cyclase",
                "AXA20105.1",
                set(),
            ]
        },
        {
            "IRqRpDzrGB9UhHJD6AzDq_6Xupj00Nte": [
                "all-trans-phytoene synthase",
                "AXA20106.1",
                set(),
            ]
        },
        {
            "ewjVum5PbpEJA4rl-BfnCAypKl5HXb7x": [
                "hypothetical protein",
                "AXA20107.1",
                set(),
            ]
        },
    ]
    assert (
        sg_object.proteins["Tdc2m3PRLsyEzjwyux6BF4arDy2mQ_Bl"].sequence
        == "MQEYCRLRRKLTLELSPEDAADVAQEAFERTLRYMRKHDGRVASPVGLLVRIALNLQIDRGRRRKHLPTALDESWEHPRWDITPEDEVVGRQSVTQLVETLDKLAPRRREAFVLCRLHGLTYQDAAKKMGIRPSVVREYLVDAVRACRDSVDWAVSRVKCNTPLVNGLELSRSG"
    )
    assert list(sg_object.assemblies.keys()) == ["lagriamide_mibig_bgc0001946"]
    assert (
        list(sg_object.assemblies["lagriamide_mibig_bgc0001946"].loci.keys())[0]
        == locus_key
    )
    assert (
        len(
            sg_object.assemblies["lagriamide_mibig_bgc0001946"].loci[locus_key].features
        )
        == 22
    )
    # check that features have their start
    start_check = [
        i.start
        for i in sg_object.assemblies["lagriamide_mibig_bgc0001946"]
        .loci[locus_key]
        .features
    ]
    start_check.sort()
    hasher = hashlib.sha256()
    hasher.update(str(start_check).encode("utf-8"))
    assert (
        hasher.hexdigest()
        == "9fbc2c80b812bcf190419da768ab4c40e4a10b1bea1326b601f0c72a86fe6cd7"
    )


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
        k: [v.description, v.other_id, v.domains]
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


def test_fasta_string_parse():
    sg_object = SocialGene()
    sg_object.parse_fasta_string(">asdads\n dasfa")
    assert list(sg_object.proteins.keys())[0] == "lZA5w9NxutVBhoqjqfsX_GX_dfugQZLQ"
    protein_parse_results = [
        {k: [v.description, v.other_id, v.domains]}
        for k, v in sg_object.proteins.items()
    ]
    assert protein_parse_results == [
        {"lZA5w9NxutVBhoqjqfsX_GX_dfugQZLQ": ["asdads", "asdads", set()]}
    ]
