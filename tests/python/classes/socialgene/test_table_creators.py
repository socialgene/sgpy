import os

from socialgene.base.socialgene import SocialGene

FIXTURE_DIR = os.path.dirname(os.path.realpath(__file__))
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)
FIXTURE_DIR = os.path.join(FIXTURE_DIR, "data", "test_genomes")
FIXTURE_DIR = os.path.join(FIXTURE_DIR, "lagriamide_mibig_bgc0001946.gbk")


def test_protein_info_table_1():
    expected_protein_info_table = [
        ("37C2BC6B34544B2E", "AXA20086.1", "sigma-70 RpoE", 174),
        ("24885C21E4F10395", "AXA20087.1", "competence protein ComEC", 64),
        ("5DB7DF559DF798B9", "AXA20088.1", "transposase", 36),
        ("380F00AA06102D54", "AXA20089.1", "hypothetical protein", 96),
        ("1CF65176FCF98848", "AXA20090.1", "hybrid trans-AT PKS/NRPS LgaA", 3553),
        ("3BFF8AA6C1E11053", "AXA20091.1", "hybrid trans-AT PKS/NRPS LgaB", 6405),
        ("5A6503F0ADBBF5DD", "AXA20092.1", "trans-AT PKS LgaC", 6799),
        ("B5363BCBA49EA647", "AXA20093.1", "trans-AT PKS LgaD", 1279),
        ("93E706B4F7EB53AF", "AXA20094.1", "acyltransferase LgaE", 370),
        ("0FB2B60341C61380", "AXA20095.1", "enoylreductase LgaF", 461),
        ("53C03058C526F0D2", "AXA20096.1", "trans-AT PKS LgaG", 8904),
        (
            "64F56D5193D6A584",
            "AXA20097.1",
            "nuclear transport factor 2 (NTF2)-like protein LgaL",
            131,
        ),
        ("A4115D03F2A9920B", "AXA20098.1", "MATE family efflux transporter LgaH", 457),
        ("631C4C8531E0AE0F", "AXA20099.1", "acylhydrolase LgaI", 337),
        ("00CD39A87B4E5579", "AXA20100.1", "transposase", 128),
        ("D9B06B39B6ED339A", "AXA20101.1", "cytochrome P450 LgaJ", 394),
        ("FF58B5837B0F03BE", "AXA20102.1", "LgaM", 244),
        ("C16214B7E088BFA5", "AXA20103.1", "ketoreductase LgaK", 273),
        (
            "2CC09AA68B15064E",
            "AXA20104.1",
            "phosphoenolpyruvate-protein phosphotransferase",
            587,
        ),
        ("F0D86567887A4777", "AXA20105.1", "squalene--hopene cyclase", 679),
        ("5CD80CF74497D50C", "AXA20106.1", "all-trans-phytoene synthase", 281),
        ("E3DBD3CB4B0B011F", "AXA20107.1", "hypothetical protein", 155),
    ]
    sg_obj = SocialGene()
    sg_obj.parse(FIXTURE_DIR)
    sg_obj.proteins
    assert [i for i in sg_obj.protein_info_table()] == expected_protein_info_table
