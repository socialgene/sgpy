import os

from socialgene.base.socialgene import SocialGene

FIXTURE_DIR = os.path.dirname(os.path.realpath(__file__))
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)

FIXTURE_DIR = os.path.join(FIXTURE_DIR, "data")

gbk_path = os.path.join(FIXTURE_DIR, "lagriamide_mibig_bgc0001946.gbk")


def test_export_locus_to_protein_1():
    sg_object = SocialGene()
    sg_object.parse(gbk_path)
    locus_key = [
        key
        for key, value in sg_object.assemblies[
            "lagriamide_mibig_bgc0001946"
        ].loci.items()
        if "bgc0001946" in key.lower()
    ][0]
    temp = list(
        sg_object.assemblies["lagriamide_mibig_bgc0001946"].loci[locus_key].features
    )
    assert sum([i.feature_is_protein() for i in temp]) == 22
    removed = sg_object.drop_non_protein_features(return_removed=True)
    assert removed == {
        "lagriamide_mibig_bgc0001946": {locus_key: {"retained": 22, "removed": 0}}
    }
    temp = list(
        sg_object.assemblies["lagriamide_mibig_bgc0001946"].loci[locus_key].features
    )
    assert len([i.feature_is_protein() for i in temp]) == 22
    assert sum([i.feature_is_protein() for i in temp]) == 22


def test_export_locus_to_protein_2():
    temp = SocialGene()
    temp.add_assembly("myassembly")
    temp.assemblies["myassembly"].add_locus(
        "my_locus",
    )
    temp.assemblies["myassembly"].loci["my_locus"].add_feature(
        protein_hash="feature_id1",
        type="protein",
        start=1,
        end=10,
        strand=1,
    )
    temp.assemblies["myassembly"].loci["my_locus"].add_feature(
        protein_hash="feature_id2",
        type="not_a_prot",
        start=1,
        end=10,
        strand=1,
    )
    assert len(temp.assemblies["myassembly"].loci["my_locus"].features) == 2
    removed = temp.drop_non_protein_features(return_removed=True)
    assert removed == {"myassembly": {"my_locus": {"retained": 1, "removed": 1}}}
    z = temp.assemblies["myassembly"].loci["my_locus"].features.pop()
    assert z.protein_hash == "feature_id1"
    removed = temp.drop_non_protein_features(return_removed=True)
    assert removed == {"myassembly": {"my_locus": {"retained": 0, "removed": 0}}}
