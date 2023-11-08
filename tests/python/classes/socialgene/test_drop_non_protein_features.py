import os

from socialgene.base.socialgene import SocialGene

FIXTURE_DIR = os.path.dirname(os.path.realpath(__file__))
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)

FIXTURE_DIR = os.path.join(FIXTURE_DIR, "data", "test_genomes")

gbk_path = os.path.join(FIXTURE_DIR, "lagriamide_mibig_bgc0001946.gbk")


def test_drop_all_non_protein_features():
    temp = SocialGene()
    temp.add_assembly("myassembly")
    temp.assemblies["myassembly"].add_locus(
        external_id="my_locus",
    )
    temp.assemblies["myassembly"].loci["my_locus"].add_feature(
        uid="feature_id1",
        type="protein",
        start=1,
        end=10,
        strand=1,
    )
    temp.assemblies["myassembly"].loci["my_locus"].add_feature(
        uid="feature_id2",
        type="not_a_prot",
        start=90,
        end=100,
        strand=-1,
    )
    assert len(temp.assemblies["myassembly"].loci["my_locus"].features) == 2
    temp.drop_all_non_protein_features()
    assert len(temp.assemblies["myassembly"].loci["my_locus"].features) == 1
    assert temp.assemblies["myassembly"].loci["my_locus"].features.pop().end == 10
