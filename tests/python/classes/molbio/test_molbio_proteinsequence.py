import pytest

from socialgene.base.molbio import ProteinSequence


def test_ProteinSequence_1():
    temp = ProteinSequence(sequence="ARNDCQEGHILKMFPSTWYVXZJU")
    assert temp.uid == "0hMjYRUCOMiDkJnVKlZ4QVMGhG8mkwdb"
    assert temp.sequence == "ARNDCQEGHILKMFPSTWYVXZJU"


def test_ProteinSequence_2():
    temp = ProteinSequence(sequence="ZJUAARNDCQEGHIILKMFPSTWTYVXZJU")
    assert (
        temp._amino_acid_count()
        == "2-1-1-1-1-1-1-1-1-2-1-1-1-1-1-1-2-1-1-1-1-2-2-2-0-0-0"
    )
    temp = ProteinSequence(sequence="STWTYVZJNDCQEGHIILKUAARMFPXZJU")
    assert (
        temp._amino_acid_count()
        == "2-1-1-1-1-1-1-1-1-2-1-1-1-1-1-1-2-1-1-1-1-2-2-2-0-0-0"
    )


def test_ProteinSequence_3():
    temp = ProteinSequence(uid="")
    assert (
        temp._amino_acid_count()
        == "0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0"
    )
    temp = ProteinSequence(sequence="")
    assert (
        temp._amino_acid_count()
        == "0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0"
    )
    temp = ProteinSequence(sequence="")
    temp.sequence = None
    assert (
        temp._amino_acid_count()
        == "0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0"
    )


def test_ProteinSequence_4():
    temp = ProteinSequence(sequence="ARNDCQEGHILKMFPSTWYVXZJU")
    assert temp.seq_len == 24


def test_fail():
    with pytest.raises(ValueError):
        ProteinSequence()


def test_fail2():
    with pytest.raises(ValueError):
        ProteinSequence(sequence="1ARNDCQEGHILKMFPSTWYVXZJU")


def test_dict():
    temp = ProteinSequence(sequence="ZJUAARNDCQEGHIILKMFPSTWTYVXZJU")
    assert temp.all_attributes() == {
        "uid": "DeAsyYJZ1XS1-LHOf15YN7OMk-gzJgbL",
        "sequence": "ZJUAARNDCQEGHIILKMFPSTWTYVXZJU",
    }


def test_crc64():
    temp = ProteinSequence(sequence="ZJUAARNDCQEGHIILKMFPSTWTYVXZJU")
    assert temp.crc64 == "29CAA4E569304732"


def test_md5():
    temp = ProteinSequence(sequence="ZJUAARNDCQEGHIILKMFPSTWTYVXZJU")
    assert temp.md5 == "e76c6fd452bb4e1376ed2f825e5d6a02"
