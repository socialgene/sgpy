import os
import tempfile
from uuid import UUID

import pytest

import socialgene.base.molbio as molbio
from socialgene.base.socialgene import SocialGene

databases = {
    "swissprot": ">sp|Q52500|THSB_THEKO Thermosome subunit beta OS=Thermococcus kodakarensis (strain ATCC BAA-918 / JCM 12380 / KOD1) (Pyrococcus kodakaraensis (strain KOD1)) OX=69014 GN=thsB PE=3 SV=1\nMAQLAGQPVVILPEGTQRYVGRDAQRLNILAARIIAETVRTTLGPKGMDKMLVDSLGDIV",
    "trembl": ">tr|Q5JCW8|Q5JCW8_THEKO Hypothetical membrane protein, conserved OS=Thermococcus kodakarensis (strain ATCC BAA-918 / JCM 12380 / KOD1) (Pyrococcus kodakaraensis (strain KOD1)) OX=69014 GN=TK0357 PE=4 SV=1\nMNPIWNL",
    "pdb": ">4D2I_1|Chains A, B|HERA|SULFOLOBUS SOLFATARICUS (273057)\nMIIGYVIGQATTQEALILA",
    "t1": "> test1\nMAQLAGQPVVILPEGTQRYVGRDAQRLNILAARIIAETVRTTLGPKGMDKMLVDSLGDIV",
    "t2": "> test 2\nMAQLAGQPVVILPEGTQRYVGRDAQRLNILAARIIAETVRTTLGPKGMDKMLVDSLGDIV",
    "t3": "> test 3 \nMAQLAGQPVVILPEGTQRYVGRDAQRLNILAARIIAETVRTTLGPKGMDKMLVDSLGDIV",
    "t4": "> test 4 \nMAQLAGQPVVILPEGTQRYVGRDAQRLNILAARIIAETVRTTLGPKGMDKMLVDSLGDIV",
    "t5": "> test 5 \nMAQLAGQPVVILPEGTQRYVGRDAQRLNILAARIIAETVRTTLGPKGMDKMLVDSLGDIV",
    "t6": "> test 6 \nmaqlagqpvvilpegtqryvgrdaqrlnilaariiaetvrttlgpkgmdkmlvdslgdiv",
    "t7": "> test|7 \nmaqlagqpvvilpegtqryvgrdaqrlnilaariiaetvrttlgpkgmdkmlvdslgdiv",
    "t8": "> te232@#4'4./'\\st|7 \nmaqlagqpvvilpegtqryvgrdaqrlnilaariiaetvrttlgpkgmdkmlvdslgdiv",
}

protein_hashes = {
    "swissprot": "r_qK-VXr_gDUb3NluH3qmgG1I09eOIOT",
    "trembl": "bhdSZvyUK6GzqHTXz4VjgtbefefmHI6x",
    "pdb": "57gtPgVt6c0hX3_B4eociVZBcgTe0TJn",
    "t1": "r_qK-VXr_gDUb3NluH3qmgG1I09eOIOT",
    "t2": "r_qK-VXr_gDUb3NluH3qmgG1I09eOIOT",
    "t3": "r_qK-VXr_gDUb3NluH3qmgG1I09eOIOT",
    "t4": "r_qK-VXr_gDUb3NluH3qmgG1I09eOIOT",
    "t5": "r_qK-VXr_gDUb3NluH3qmgG1I09eOIOT",
    "t6": "r_qK-VXr_gDUb3NluH3qmgG1I09eOIOT",
    "t7": "r_qK-VXr_gDUb3NluH3qmgG1I09eOIOT",
    "t8": "r_qK-VXr_gDUb3NluH3qmgG1I09eOIOT",
}


@pytest.mark.parametrize(
    "k,v", [(k, v) for k, v in databases.items()], ids=[k for k, v in databases.items()]
)
def test_assembly_from_file(k, v):
    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as f:
        f.write(v)
        f.seek(0)
        a = SocialGene()
        a.parse(f.name)
        list(a.assemblies.keys())[0] == os.path.basename(f.name)


@pytest.mark.parametrize(
    "k,v", [(k, v) for k, v in databases.items()], ids=[k for k, v in databases.items()]
)
def test_assembly_from_string(k, v):
    a = SocialGene()
    a.parse(v)
    # test that the ids are a uuid4
    assert UUID(list(a.assemblies.keys())[0], version=4)


@pytest.mark.parametrize(
    "k,v", [(k, v) for k, v in databases.items()], ids=[k for k, v in databases.items()]
)
def test_protein_hash_from_file(k, v):
    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as f:
        f.write(v)
        f.seek(0)
        a = SocialGene()
        a.parse(f.name)
        # test that the ids are a uuid4
        assert list(a.proteins.keys())[0] == protein_hashes[k]


@pytest.mark.parametrize(
    "k,v", [(k, v) for k, v in databases.items()], ids=[k for k, v in databases.items()]
)
def test_protein_hash_from_string(k, v):
    a = SocialGene()
    a.parse(v)
    assert list(a.proteins.keys())[0] == protein_hashes[k]


def test_parse_multiple():
    a = SocialGene()
    for k, v in databases.items():
        a.parse(v)
    len(a.proteins) == 3


def test_parse_defline_magic_true():
    a = SocialGene()
    a.parse("\n".join(databases.values()), defline_magic=True)
    len(a.proteins) == 3
    assert {k: v.all_attributes() for k, v in a.proteins.items()} == {
        "r_qK-VXr_gDUb3NluH3qmgG1I09eOIOT": {
            "description": "sp|Q52500|THSB_THEKO Thermosome subunit beta OS=Thermococcus kodakarensis (strain ATCC BAA-918 / JCM 12380 / KOD1) (Pyrococcus kodakaraensis (strain KOD1)) OX=69014 GN=thsB PE=3 SV=1",
            "domains": set(),
            "external_id": "Q52500",
            "seq_len": 60,
        },
        "bhdSZvyUK6GzqHTXz4VjgtbefefmHI6x": {
            "description": "tr|Q5JCW8|Q5JCW8_THEKO Hypothetical membrane protein, conserved OS=Thermococcus kodakarensis (strain ATCC BAA-918 / JCM 12380 / KOD1) (Pyrococcus kodakaraensis (strain KOD1)) OX=69014 GN=TK0357 PE=4 SV=1",
            "domains": set(),
            "external_id": "Q5JCW8",
            "seq_len": 7,
        },
        "57gtPgVt6c0hX3_B4eociVZBcgTe0TJn": {
            "description": "4D2I_1|Chains A, B|HERA|SULFOLOBUS SOLFATARICUS (273057)",
            "domains": set(),
            "external_id": "4D2I_1|Chains",
            "seq_len": 19,
        },
    }


def test_parse_defline_magic_false():
    a = SocialGene()
    a.parse("\n".join(databases.values()), defline_magic=False)
    len(a.proteins) == 3
    assert {k: v.all_attributes() for k, v in a.proteins.items()} == {
        "r_qK-VXr_gDUb3NluH3qmgG1I09eOIOT": {
            "description": "sp|Q52500|THSB_THEKO Thermosome subunit beta OS=Thermococcus kodakarensis (strain ATCC BAA-918 / JCM 12380 / KOD1) (Pyrococcus kodakaraensis (strain KOD1)) OX=69014 GN=thsB PE=3 SV=1",
            "domains": set(),
            "external_id": "sp|Q52500|THSB_THEKO",
            "seq_len": 60,
        },
        "bhdSZvyUK6GzqHTXz4VjgtbefefmHI6x": {
            "description": "tr|Q5JCW8|Q5JCW8_THEKO Hypothetical membrane protein, conserved OS=Thermococcus kodakarensis (strain ATCC BAA-918 / JCM 12380 / KOD1) (Pyrococcus kodakaraensis (strain KOD1)) OX=69014 GN=TK0357 PE=4 SV=1",
            "domains": set(),
            "external_id": "tr|Q5JCW8|Q5JCW8_THEKO",
            "seq_len": 7,
        },
        "57gtPgVt6c0hX3_B4eociVZBcgTe0TJn": {
            "description": "4D2I_1|Chains A, B|HERA|SULFOLOBUS SOLFATARICUS (273057)",
            "domains": set(),
            "external_id": "4D2I_1|Chains",
            "seq_len": 19,
        },
    }


def test_parse_2():
    a = SocialGene()
    a.parse(
        "\n".join([v for k, v in databases.items() if k in ["swissprot", "trembl"]]),
        defline_magic=True,
    )
    for ak, av in a.assemblies.items():
        for lk, lv in av.loci.items():
            for fk in lv.features:
                assert isinstance(fk.parent_object, molbio.Locus)
                del fk.parent_object
                if fk.external_id == "Q5JCW8":
                    assert fk.all_attributes() == {
                        "description": "tr|Q5JCW8|Q5JCW8_THEKO Hypothetical membrane protein, conserved OS=Thermococcus kodakarensis (strain ATCC BAA-918 / JCM 12380 / KOD1) (Pyrococcus kodakaraensis (strain KOD1)) OX=69014 GN=TK0357 PE=4 SV=1",
                        "external_id": "Q5JCW8",
                        "frameshifted": None,
                        "goterms": None,
                        "incomplete": None,
                        "internal_stop": None,
                        "locus_tag": None,
                        "missing_C_terminus": None,
                        "missing_N_terminus": None,
                        "missing_start": None,
                        "missing_stop": None,
                        "note": None,
                        "partial_in_the_middle_of_a_contig": None,
                        "partial_on_complete_genome": None,
                        "too_short_partial_abutting_assembly_gap": None,
                        "type": "protein",
                        "uid": "bhdSZvyUK6GzqHTXz4VjgtbefefmHI6x",
                    }
                if fk.external_id == "Q52500":
                    assert fk.all_attributes() == {
                        "description": "sp|Q52500|THSB_THEKO Thermosome subunit beta OS=Thermococcus kodakarensis (strain ATCC BAA-918 / JCM 12380 / KOD1) (Pyrococcus kodakaraensis (strain KOD1)) OX=69014 GN=thsB PE=3 SV=1",
                        "external_id": "Q52500",
                        "frameshifted": None,
                        "goterms": None,
                        "incomplete": None,
                        "internal_stop": None,
                        "locus_tag": None,
                        "missing_C_terminus": None,
                        "missing_N_terminus": None,
                        "missing_start": None,
                        "missing_stop": None,
                        "note": None,
                        "partial_in_the_middle_of_a_contig": None,
                        "partial_on_complete_genome": None,
                        "too_short_partial_abutting_assembly_gap": None,
                        "type": "protein",
                        "uid": "r_qK-VXr_gDUb3NluH3qmgG1I09eOIOT",
                    }
