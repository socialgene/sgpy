import pytest

from socialgene.base.molbio import ProteinSequence


def test_ProteinSequence_1():
    temp = ProteinSequence(sequence="ARNDCQEGHILKMFPSTWYVXZJU")
    assert temp.hash_id == "0hMjYRUCOMiDkJnVKlZ4QVMGhG8mkwdb"
    assert temp.sequence == "ARNDCQEGHILKMFPSTWYVXZJU"
    temp.hash()
    assert temp.hash_id == "0hMjYRUCOMiDkJnVKlZ4QVMGhG8mkwdb"


def test_ProteinSequence_2():
    temp = ProteinSequence(sequence="ZJUAARNDCQEGHIILKMFPSTWTYVXZJU")
    assert (
        temp.count_amino_acids()
        == "2-1-1-1-1-1-1-1-1-2-1-1-1-1-1-1-2-1-1-1-1-2-2-2-0-0-0"
    )


def test_ProteinSequence_3():
    temp = ProteinSequence(sequence="ARNDCQEGHILKMFPSTWYVXZJU")
    assert temp.sequence_length() == 24


def test_fail():
    with pytest.raises(ValueError):
        ProteinSequence()


def test_fail2():
    with pytest.raises(ValueError):
        ProteinSequence(sequence="1ARNDCQEGHILKMFPSTWYVXZJU")
