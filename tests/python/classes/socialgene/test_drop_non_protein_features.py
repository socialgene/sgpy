import os
from socialgene.base.socialgene import SocialGene


FIXTURE_DIR = os.path.dirname(os.path.realpath(__file__))
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)

FIXTURE_DIR = os.path.join(FIXTURE_DIR, "data")

gbk_path = os.path.join(FIXTURE_DIR, "lagriamide_mibig_bgc0001946.gbk")


def test_export_locus_to_protein():
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
