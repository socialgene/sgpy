from unittest.mock import patch

from socialgene.base.molbio import Assembly, Molbio


def test_empty():
    assert Molbio().__dict__ == {"assemblies": {}, "proteins": {}}


def test_empty2():
    assert Molbio().get_all_feature_uids() == []


def add_assembly():
    a = Molbio()
    a.add_assembly(uid="id1")
    assert isinstance(a.assemblies["id1"], Assembly)
    assert a.assemblies["id1"].name == "id1"


@patch("socialgene.base.molbio.uuid4")
def test_add_assembly_uuid(my_method):
    my_method.return_value = "ASDSDSD"
    a = Molbio()
    a.add_assembly()
    assert list(a.assemblies.keys())[0] == "ASDSDSD"
