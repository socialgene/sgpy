import os
import tempfile
from pathlib import Path

import pytest

from socialgene.utils.goterms import (
    AltRel,
    GoNode,
    GoRelationship,
    extract_goterm_int_as_str,
    parse,
)

FIXTURE_DIR = os.path.dirname(os.path.realpath(__file__))
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)
gofile = Path(FIXTURE_DIR, "data", "obo", "go-basic-head.obo")


def test_extract_goterm_int_as_str():
    assert extract_goterm_int_as_str("GO:1904659") == "1904659"
    # pulling from an unknown longer string so it truncates, no great way to kno termination
    assert extract_goterm_int_as_str("GO:19046590") == "1904659"
    assert extract_goterm_int_as_str("GO:190465") is None
    assert extract_goterm_int_as_str("1904659") is None
    assert extract_goterm_int_as_str("") is None


def test_GoNode_slots():
    assert GoNode().__slots__ == ["uid", "name", "namespace", "definition"]


def test_GoNode_add_id():
    with pytest.raises(Exception):
        a = GoNode()
        a.add_id("")
    with pytest.raises(Exception):
        a = GoNode()
        a.add_id("GO:190465")
    a = GoNode()
    a.add_id("GO:1904659")
    assert a.uid == "1904659"


def test_GoNode_add_name():
    a = GoNode()
    a.add_name("test1")
    assert a.name == "test1"
    a.add_name("name:test2")
    assert a.name == "test2"
    a.add_name("name: test3 ")
    assert a.name == "test3"


def test_GoNode_add_namespace():
    a = GoNode()
    a.add_namespace("test1")
    assert a.namespace == "test1"
    a.add_namespace("namespace:test2")
    assert a.namespace == "test2"
    a.add_namespace("namespace: test3 ")
    assert a.namespace == "test3"


def test_GoNode_add_definition():
    a = GoNode()
    a.add_definition("test1")
    assert a.definition == "test1"
    a.add_definition("definition:test2")
    assert a.definition == "test2"
    a.add_definition("definition: test3 ")
    assert a.definition == "test3"


def test_assign_output():
    a = GoNode()
    for i in ["id:GO:1904659", "name: 2", "namespace: 3", "definition:4"]:
        a.assign(i)
    assert a.output() == ("1904659", "2", "3", "4")


def test_GoRelationship_add_end():
    a = GoRelationship()
    a.add_end("GO:1904659")
    assert a.end == "1904659"


def test_GoRelationship_add_type():
    a = GoRelationship()
    a.add_type("relationship: regulates GO:0006310 ! DNA recombination")
    assert a.type == "REGULATES"
    a.add_type("relationship: ahh!!! GO:0006310 ! DNA recombination")
    assert a.type == "AHH!!!"
    a.add_type("is_a: GO:0051052 ! regulation of DNA metabolic process")
    assert a.type == "IS_A"


def test_GoRelationship_assign():
    a = GoRelationship()
    a.assign("2", "relationship: regulates GO:0006310 ! DNA recombination")
    assert a.type == "REGULATES"
    assert a.start == "2"
    assert a.output() == ("2", "0006310", "REGULATES")
    a.assign("1", "is_a: GO:0051052 ! regulation of DNA metabolic process")
    assert a.start == "1"
    assert a.type == "IS_A"
    assert a.output() == ("1", "0051052", "IS_A")


def test_AltRel():
    a = AltRel()
    a.assign(main_id="1", line="alt_id: GO:0000013")
    assert a.output() == ("1", "0000013", "ALTERNATE")


def test_parse():
    with tempfile.TemporaryDirectory() as tmpdirname:
        parse(
            outdir=tmpdirname,
            filepath=gofile,
        )
        with open(Path(tmpdirname, "goterms"), "r") as h:
            assert h.readlines() == [
                "0000001\tmitochondrion inheritance\tbiological_process\t\n",
                "0000002\tmitochondrial genome maintenance\tbiological_process\t\n",
                "0000016\tlactase activity\tbiological_process\t\n",
                "0000018\tregulation of DNA recombination\tbiological_process\t\n",
                "0000013\tobsolete thioredoxin\tmolecular_function\t\n",
                "0000008\tobsolete thioredoxin\tmolecular_function\t\n",
            ]

        with open(Path(tmpdirname, "goterm_edgelist"), "r") as h:
            assert h.readlines() == [
                "0000001\t0048311\tIS_A\n",
                "0000002\t0007005\tIS_A\n",
                "0000016\t0042946\tIS_A\n",
                "0000018\t0006310\tREGULATES\n",
                "0000008\t0000013\tALTERNATE\n",
            ]


def test_parse_fail_args():
    with pytest.raises(ValueError):
        parse(outdir="", filepath=None, url=None)
