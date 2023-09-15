import os
import tempfile
from collections import OrderedDict
from pathlib import Path

import pytest

from socialgene.hashing.hashing import hasher
from socialgene.parsers.hmmmodel import (
    HMM_SOURCES,
    HmmModel,
    HmmParse,
    parse_hmmlist_input,
)

tmpdir = tempfile.TemporaryDirectory()


FIXTURE_DIR = os.path.dirname(os.path.realpath(__file__))
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)
FIXTURE_DIR = Path(FIXTURE_DIR, "data", "test_genomes")
gbk_path = Path(FIXTURE_DIR, "lagriamide_mibig_bgc0001946.gbk")


FIXTURE_DIR = os.path.dirname(os.path.realpath(__file__))
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)
hmm_dir = Path(FIXTURE_DIR, "data", "hmms")


def test_HmmModel():
    a = HmmModel()
    assert a.__dict__ == OrderedDict(
        [
            ("id", "None_None"),
            ("source", None),
            ("rel_path", None),
            ("name", None),
            ("acc", None),
            ("notes", ""),
            ("description", None),
            ("date", None),
            ("hash", None),
            ("hash_used", None),
            ("model_length", 0),
            ("category", None),
            ("subcategory", None),
            ("ga", None),
            ("tc", None),
            ("nc", None),
        ]
    )


def test_notes_1():
    a = HmmModel()
    a._notes = ["1", "a"]
    assert a.__dict__ == OrderedDict(
        [
            ("id", "None_None"),
            ("source", None),
            ("rel_path", None),
            ("name", None),
            ("acc", None),
            ("notes", "1; a"),
            ("description", None),
            ("date", None),
            ("hash", None),
            ("hash_used", None),
            ("model_length", 0),
            ("category", None),
            ("subcategory", None),
            ("ga", None),
            ("tc", None),
            ("nc", None),
        ]
    )


def test_notes_2():
    a = HmmModel()
    a._notes = "1"
    assert a.__dict__ == OrderedDict(
        [
            ("id", "None_None"),
            ("source", None),
            ("rel_path", None),
            ("name", None),
            ("acc", None),
            ("notes", "1"),
            ("description", None),
            ("date", None),
            ("hash", None),
            ("hash_used", None),
            ("model_length", 0),
            ("category", None),
            ("subcategory", None),
            ("ga", None),
            ("tc", None),
            ("nc", None),
        ]
    )


def test_parse_hmmlist_input():
    assert parse_hmmlist_input("a") == []


@pytest.mark.parametrize("n", HMM_SOURCES)
def test_parse_hmmlist_individual_input(n):
    assert parse_hmmlist_input(n) == [n]


def test_parse_hmmlist_all_input():
    assert parse_hmmlist_input("all") == HMM_SOURCES


def test_model_strings():
    a = HmmParse()
    hmm_path = Path(hmm_dir, "pks.hmm")
    a.read(hmm_path, hmm_path.parent)
    assert {k: hasher(v.model_string()) for k, v in a.models.items()} == {
        0: "UVT0GljnP1J2v8JCa8NQvgbOdNACNAMW",
        1: "LfWmuMMJjYIqAV3m97uBS8THZvQostqB",
        2: "WFr9xOaMtY0FBJsbyVT8JbGG20QS-3Mz",
        3: "7MlcIMFQlyAPPM_3v8JPiGReEUZtmH-j",
        4: "tgDdTsE0JwtXje03tmScunq6u-90Ddsa",
        5: "NuKaOPo-E2CA_Y3R0oihDGFtEyM-IvOR",
        6: "V3XYonn5bh9i1cAMxjTqDQjncnCYWc3S",
        7: "uQOY4LWrFFvdYvwWyubenokylzZ1nMIx",
        8: "FnGLSiT418oxwd8gQzOGj36FEEdOpkoV",
        9: "0rC83aj5OuHkGcFJ6zSf9IhseGp2ui4c",
        10: "VSUeZMI6jU1IsMmzHM83r22fRLAwxyQC",
        11: "16NLuotpFujFlJI0E1DxOD2zX71CDhFq",
        12: "i2lvejQGdcV7tAmCLNykHnFas9mx4VOZ",
        13: "bh1uOdotjJffuzTsM7y9E5qsBHZR49WI",
        14: "S6o9dIJ7SapVtG249R1m2S7X5lQFoJMw",
        15: "DrXz7D4smiGK6WgO9FmBZZldKdxqxNmv",
        16: "e7aomAooAwKxwexbr8pakaSPLhFlafnD",
        17: "E63qgM2y9ny8xRtuDoDcnVSrJ6pTP4qu",
        18: "WNnUzK7i3T-XhsSYKiABflSeW4jJ4pCd",
        19: "Rns0jpXjgNdyU5o_7Dg5bA43ZN1acVEO",
        20: "HgZNGEsy3y3TpPYRyUXsp92A1uLvYGMc",
        21: "6aikW43Annzy4BBFOUE7QBGYN852TznL",
        22: "hadncr1mUju9K3WiMP0qNGsVHKYRhnsE",
        23: "Y12LiFgPg7IN4DhB1ESenlWSv3DZz2Yc",
        24: "Xl-4khMnF147PbxoAry539dJTF2v36Hw",
        25: "Zjgs6wCgJBfsiWwCGEi9JSXi_sD0Y-jb",
    }


def test_model_strings_hash_as_name():
    a = HmmParse()
    hmm_path = Path(hmm_dir, "pks.hmm")
    a.read(hmm_path, hmm_path.parent)
    assert {
        k: hasher(v.model_string(hash_as_name=True)) for k, v in a.models.items()
    } == {
        0: "nbuSy-HMEOdfqNEPS9QvUBgkKTXI_DT0",
        1: "3sku308mR2U7KMQK3S9y8FGJLzR3B6P7",
        2: "gtOAw7gikXfLJt5x_MbaEnJua_Be7XUI",
        3: "9ufSNpXeD__WdrSyn_x9I6j-OC1aYr4N",
        4: "RmOqZsyilo5L8fDV6ktcrgPsIU-_A-4q",
        5: "2oW3DCN2lRSGlC4QkU5zdX54Lmgdb5ZK",
        6: "3CJJ7GrDal-6U55-T8_00xRFsm647G60",
        7: "qIbEq0YMNC2hK4KH5uuNam9crNhwkARz",
        8: "PrnJIUd3XzsfkcF6DCl6URbEHiZ1GaVO",
        9: "SmpT33h7e5YNDjDxbc14DwHtpmn1fxqw",
        10: "BpDePy5sT8NqS-sAWb4PzgArI6ea19U_",
        11: "TMx0W9jzUdGC_5ObvijNOoCXwVbMs7Xw",
        12: "Wt6ezgaDACWLbZWPoKsF2rkPzY57EzGY",
        13: "PzJ426PpyYMCyUqQBrwL2c-Hgr_Q2cS6",
        14: "6FGiUyghLjH1pTXnSSIxCkHm8IDthb0t",
        15: "JBzkrQGZK9lpZtGpaQzWvOgjNqRrvFM9",
        16: "fdLmvUyQB5kMZBWlXwDPiRYc2dJ5BcQs",
        17: "Cq-JeXOxFIp2lLfhVLWI96WVEoWoN61K",
        18: "nwJMlrNj5BOWV3ofYEOU97lRo3osyfBP",
        19: "yfBLCeh4Yg6Iy2QNBRcPKTARJTSocfLf",
        20: "oFLjAvhouQEraFhld9-Hbiv0dLl5ObvV",
        21: "GwFrMc0kSa56HCnW0swcrubq8fws53Sq",
        22: "mW1PWwrBR0L-IPoPSZyrPB1ePbzgP1fq",
        23: "yAF8kBz2Qa2TDXEMjjImolhYMHAfBLgj",
        24: "iyI7QzrmVMhSqz0QsuECQ2IOQBS3yYeW",
        25: "toztathlMpp3XOqhch-9isoJ0DHQG7IW",
    }


def test_notes_3():
    a = HmmModel()
    a._notes = "1"
    assert a.__dict__ == OrderedDict(
        [
            ("id", "None_None"),
            ("source", None),
            ("rel_path", None),
            ("name", None),
            ("acc", None),
            ("notes", "1"),
            ("description", None),
            ("date", None),
            ("hash", None),
            ("hash_used", None),
            ("model_length", 0),
            ("category", None),
            ("subcategory", None),
            ("ga", None),
            ("tc", None),
            ("nc", None),
        ]
    )


def test_model_write():
    a = HmmParse()
    hmm_path = Path(hmm_dir, "pks.hmm")
    a.read(hmm_path, hmm_path.parent)
    with tempfile.NamedTemporaryFile() as fp:
        for i in a.models.values():
            i.write(outpath=fp.name, mode="a", hash_as_name=False)
        with open(fp.name, "r") as target:
            with open(hmm_path, "r") as expected:
                a1 = target.readlines()
                a2 = expected.readlines()
                for i1, i2 in zip(a1, a2):
                    assert i1.strip() == i2.strip()
