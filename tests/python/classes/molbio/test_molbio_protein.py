import pytest
from socialgene.base.molbio import Protein


temp = Protein(
    sequence="ARNDCQEGHILKMFPSTWYVXZJU",
    description="description",
    other_id="other_id",
    domains=[],
)


def test_protein():
    assert temp.description == "description"
    assert temp.domains == []
    assert temp.hash_id == "0hMjYRUCOMiDkJnVKlZ4QVMGhG8mkwdb"
    assert temp.other_id == "other_id"
    assert temp.sequence == "ARNDCQEGHILKMFPSTWYVXZJU"


def test_fail():
    with pytest.raises(ValueError):
        Protein()
