import os
import tempfile

import socialgene.base.molbio as molbio
from socialgene.base.socialgene import SocialGene

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
        {k: [v.description, v.external_protein_id, v.domains]}
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
    # Note this wouldn't be good outside her because diff features can same id value
    # This (especially {"start":i.start, "end": i.end, "strand":i.strand, "type": i.type}) is verbose
    # so it's hopefully easier to pinpoint any arising errors
    feature_check = {
        i.protein_hash: {
            "start": i.start,
            "end": i.end,
            "strand": i.strand,
            "type": i.type,
        }
        for i in sg_object.assemblies["lagriamide_mibig_bgc0001946"]
        .loci[locus_key]
        .features
    }
    assert feature_check == {
        "9stFB1fGjCdZVZWHLVI3OD4A_DV3WcV6": {
            "start": 86717,
            "end": 87112,
            "strand": 1,
            "type": "CDS",
        },
        "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L": {
            "start": 92859,
            "end": 93680,
            "strand": 1,
            "type": "CDS",
        },
        "RyDIaUZc_b21_kQalx7J3yNO4l5f-439": {
            "start": 53400,
            "end": 57239,
            "strand": 1,
            "type": "CDS",
        },
        "Tdc2m3PRLsyEzjwyux6BF4arDy2mQ_Bl": {
            "start": 1190,
            "end": 1714,
            "strand": -1,
            "type": "CDS",
        },
        "Qi23auOUTcBzWTmDHuGinrzIuqH7-zVn": {
            "start": 88580,
            "end": 89593,
            "strand": 1,
            "type": "CDS",
        },
        "WbViYzQw8y-XfCQMgQXkedGduNMJPa14": {
            "start": 93904,
            "end": 95667,
            "strand": -1,
            "type": "CDS",
        },
        "4_8182J88axMDpFJBZI6kLNJAu8Ittm3": {
            "start": 95865,
            "end": 97904,
            "strand": 1,
            "type": "CDS",
        },
        "ptq1NGhBcUp3TIEqvAUxnnp4LOKwINvn": {
            "start": 1915,
            "end": 2109,
            "strand": 1,
            "type": "CDS",
        },
        "nk3UJUyLaWr8LochtohJ9L6eugdChZL9": {
            "start": 90318,
            "end": 90704,
            "strand": 1,
            "type": "CDS",
        },
        "IRqRpDzrGB9UhHJD6AzDq_6Xupj00Nte": {
            "start": 98026,
            "end": 98871,
            "strand": 1,
            "type": "CDS",
        },
        "-l7xLyFZbiZENPLq_GML8JyTRF1Srawr": {
            "start": 2205,
            "end": 2315,
            "strand": 1,
            "type": "CDS",
        },
        "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5": {
            "start": 33004,
            "end": 53403,
            "strand": 1,
            "type": "CDS",
        },
        "DTee9G4M8sEfnM4HaPfI37rT74pq7M_G": {
            "start": 57289,
            "end": 58401,
            "strand": 1,
            "type": "CDS",
        },
        "ewjVum5PbpEJA4rl-BfnCAypKl5HXb7x": {
            "start": 98927,
            "end": 99394,
            "strand": 1,
            "type": "CDS",
        },
        "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh": {
            "start": 59909,
            "end": 86623,
            "strand": 1,
            "type": "CDS",
        },
        "du1Ncfm5UYiFYgDWD8KW1AQJNHlAcVXL": {
            "start": 90836,
            "end": 92020,
            "strand": 1,
            "type": "CDS",
        },
        "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl": {
            "start": 3132,
            "end": 13793,
            "strand": 1,
            "type": "CDS",
        },
        "mB22-i4RqtslyO7_HappM4rJ4Z2Qbkfn": {
            "start": 58398,
            "end": 59783,
            "strand": 1,
            "type": "CDS",
        },
        "5IYMhENey2WCMrKPUz3AqBIZuFSv6DPP": {
            "start": 92083,
            "end": 92817,
            "strand": 1,
            "type": "CDS",
        },
        "MSHRSCZfdBJP8vdJdaXfeZrThH_4EUMm": {
            "start": 87147,
            "end": 88520,
            "strand": 1,
            "type": "CDS",
        },
        "T_DzOorDp3ROhRRBtuXP3xyAPorpTVD0": {
            "start": 2828,
            "end": 3118,
            "strand": 1,
            "type": "CDS",
        },
        "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ": {
            "start": 13790,
            "end": 33007,
            "strand": 1,
            "type": "CDS",
        },
    }


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
