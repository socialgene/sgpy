import socialgene
from socialgene.base.molbio import Assembly
from socialgene.base.socialgene import SocialGene

temp = SocialGene()


temp.add_protein(
    sequence="ARNDCQEGHILKMFPSTWYVXZJU",
    description="description",
    external_protein_id="external_protein_id",
    domains=[],
)


temp.add_assembly("myassembly")
temp.assemblies["myassembly"].add_locus(
    "my_locus",
)

temp.assemblies["myassembly"].loci["my_locus"].add_feature(
    protein_hash="feature_id",
    type="prot",
    start=1,
    end=10,
    strand=1,
)


def test_add_protein():
    assert temp.proteins["20F58F6F237F111D"].description == "description"
    assert temp.proteins["20F58F6F237F111D"].domains == []
    assert temp.proteins["20F58F6F237F111D"].hash_id == "20F58F6F237F111D"
    assert (
        temp.proteins["20F58F6F237F111D"].external_protein_id == "external_protein_id"
    )
    assert temp.proteins["20F58F6F237F111D"].sequence == "ARNDCQEGHILKMFPSTWYVXZJU"


def test_add_assembly():
    assert isinstance(temp.assemblies, dict)
    assert isinstance(temp.assemblies["myassembly"], socialgene.base.molbio.Assembly)
    assert list(temp.assemblies.keys()) == ["myassembly"]
    assert isinstance(list(temp.assemblies.values())[0], Assembly)
    assert len(list(temp.assemblies.values())) == 1


def test_add_locus1():
    assert isinstance(
        temp.assemblies["myassembly"].loci["my_locus"], socialgene.base.molbio.Locus
    )
    assert list(temp.assemblies["myassembly"].loci.keys()) == ["my_locus"]
    assert len(temp.assemblies["myassembly"].loci["my_locus"].features) == 1


def test_add_feature():
    assert isinstance(temp.assemblies["myassembly"].loci["my_locus"].features, set)
    a = temp.assemblies["myassembly"].loci["my_locus"].features.pop()
    assert isinstance(a, socialgene.base.molbio.Feature)
    assert a.start == 1
    assert a.end == 10
    assert a.strand == 1
    assert a.protein_hash == "feature_id"
    assert a.type == "prot"
