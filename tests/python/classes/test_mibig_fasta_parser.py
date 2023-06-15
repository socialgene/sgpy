import os
import tempfile

from socialgene.base.socialgene import SocialGene

from .test_mibig_gbk_parser import PROTEIN_DICT

DIRECTORY_OF_THIS_FILE = os.path.dirname(os.path.realpath(__file__))
FIXTURE_DIR = os.path.dirname(DIRECTORY_OF_THIS_FILE)

FIXTURE_DIR = os.path.join(FIXTURE_DIR, "data", "test_genomes")


def test_fasta_file_parse():
    # just create the fasta from the genbank file
    with tempfile.NamedTemporaryFile(suffix=".faa") as fp:
        sg_object = SocialGene()
        gbk_path = os.path.join(FIXTURE_DIR, "lagriamide_mibig_bgc0001946.gbk")
        sg_object.parse(gbk_path)
        sg_object.write_fasta(outpath=fp.name)
        fasta_object = SocialGene()
        fasta_object.parse(fp.name)

    protein_parse_results = {
        k: [v.description, v.external_protein_id, v.domains]
        for k, v in fasta_object.proteins.items()
    }
    # fmt: off
    assert protein_parse_results == {'37C2BC6B34544B2E': ['37C2BC6B34544B2E', '37C2BC6B34544B2E', set()], '24885C21E4F10395': ['24885C21E4F10395', '24885C21E4F10395', set()], '5DB7DF559DF798B9': ['5DB7DF559DF798B9', '5DB7DF559DF798B9', set()], '380F00AA06102D54': ['380F00AA06102D54', '380F00AA06102D54', set()], '1CF65176FCF98848': ['1CF65176FCF98848', '1CF65176FCF98848', set()], '3BFF8AA6C1E11053': ['3BFF8AA6C1E11053', '3BFF8AA6C1E11053', set()], '5A6503F0ADBBF5DD': ['5A6503F0ADBBF5DD', '5A6503F0ADBBF5DD', set()], 'B5363BCBA49EA647': ['B5363BCBA49EA647', 'B5363BCBA49EA647', set()], '93E706B4F7EB53AF': ['93E706B4F7EB53AF', '93E706B4F7EB53AF', set()], '0FB2B60341C61380': ['0FB2B60341C61380', '0FB2B60341C61380', set()], '53C03058C526F0D2': ['53C03058C526F0D2', '53C03058C526F0D2', set()], '64F56D5193D6A584': ['64F56D5193D6A584', '64F56D5193D6A584', set()], 'A4115D03F2A9920B': ['A4115D03F2A9920B', 'A4115D03F2A9920B', set()], '631C4C8531E0AE0F': ['631C4C8531E0AE0F', '631C4C8531E0AE0F', set()], '00CD39A87B4E5579': ['00CD39A87B4E5579', '00CD39A87B4E5579', set()], 'D9B06B39B6ED339A': ['D9B06B39B6ED339A', 'D9B06B39B6ED339A', set()], 'FF58B5837B0F03BE': ['FF58B5837B0F03BE', 'FF58B5837B0F03BE', set()], 'C16214B7E088BFA5': ['C16214B7E088BFA5', 'C16214B7E088BFA5', set()], '2CC09AA68B15064E': ['2CC09AA68B15064E', '2CC09AA68B15064E', set()], 'F0D86567887A4777': ['F0D86567887A4777', 'F0D86567887A4777', set()], '5CD80CF74497D50C': ['5CD80CF74497D50C', '5CD80CF74497D50C', set()], 'E3DBD3CB4B0B011F': ['E3DBD3CB4B0B011F', 'E3DBD3CB4B0B011F', set()]}
    assert {k: v.sequence for k, v in sg_object.proteins.items()} == PROTEIN_DICT


def test_fasta_string_parse():
    sg_object = SocialGene()
    sg_object.parse_fasta_string(">asdads\n dasfa")
    assert list(sg_object.proteins.keys())[0] == "6DD9D5BDDAC00000"
    protein_parse_results = [
        {k: [v.description, v.external_protein_id, v.domains]}
        for k, v in sg_object.proteins.items()
    ]
    assert protein_parse_results == [{"6DD9D5BDDAC00000": ["asdads", "asdads", set()]}]
