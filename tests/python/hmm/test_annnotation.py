import os

from socialgene.base.socialgene import SocialGene
from socialgene.hmm.hmmer import HMMER

FIXTURE_DIR = os.path.dirname(os.path.realpath(__file__))
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)
FIXTURE_DIR = os.path.join(FIXTURE_DIR, "data", "test_genomes")
gbk_path = os.path.join(FIXTURE_DIR, "lagriamide_mibig_bgc0001946.gbk")


FIXTURE_DIR = os.path.dirname(os.path.realpath(__file__))
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)
FIXTURE_DIR = os.path.join(FIXTURE_DIR, "data")
hmm_path = os.path.join(FIXTURE_DIR, "pks.hmm")


expected_proteins = [
    "00CD39A87B4E5579",
    "0FB2B60341C61380",
    "1CF65176FCF98848",
    "24885C21E4F10395",
    "2CC09AA68B15064E",
    "37C2BC6B34544B2E",
    "380F00AA06102D54",
    "3BFF8AA6C1E11053",
    "53C03058C526F0D2",
    "5A6503F0ADBBF5DD",
    "5CD80CF74497D50C",
    "5DB7DF559DF798B9",
    "631C4C8531E0AE0F",
    "64F56D5193D6A584",
    "93E706B4F7EB53AF",
    "A4115D03F2A9920B",
    "B5363BCBA49EA647",
    "C16214B7E088BFA5",
    "D9B06B39B6ED339A",
    "E3DBD3CB4B0B011F",
    "F0D86567887A4777",
    "FF58B5837B0F03BE",
]


def test_hmmscan():
    # have to set env here because other tests change it
    from socialgene.config import env_vars

    env_vars["HASHING_ALGORITHM"] = "crc64"

    env_vars["HMMSEARCH_IEVALUE"] = 0.1
    env_vars["HMMSEARCH_Z"] = 4
    sg_object = SocialGene()
    sg_object.parse(gbk_path)
    protein_id_list = list(sg_object.proteins.keys())

    h = HMMER()
    h.hmmpress(hmm_path)
    sg_object.annotate_with_hmmscan(
        protein_id_list=protein_id_list, hmm_filepath=hmm_path, cpus=1
    )
    prots = list(sg_object.proteins.keys())
    prots.sort()
    assert prots == expected_proteins
    # note HMM model file was created using sha512t24u which is why hashes are that
    assert sorted(
        [i.hmm_id for i in sg_object.proteins["37C2BC6B34544B2E"].domains]
    ) == ["DIJQMpAiLKGDPgcpc1IuBzFdf7FhTYu5", "xJwofaGb0EIZrSxSeZL5xS6thEM7ck7U"]
    for i in sg_object.proteins["37C2BC6B34544B2E"].domains:
        if i.hmm_id == "DIJQMpAiLKGDPgcpc1IuBzFdf7FhTYu5":
            # '{s: getattr(i, s, None) for s in i.__slots__}' turns the Domain object into a dict
            assert {s: getattr(i, s, None) for s in i.__slots__} == {
                "hmm_id": "DIJQMpAiLKGDPgcpc1IuBzFdf7FhTYu5",
                "env_from": 4,
                "env_to": 66,
                "seq_pro_score": 28.6,
                "evalue": -3,
                "i_evalue": -3,
                "domain_bias": 0.1,
                "domain_score": 26.9,
                "seq_pro_bias": 0.2,
                "hmm_from": 19,
                "hmm_to": 69,
                "ali_from": 16,
                "ali_to": 65,
                "exponentialized": True,
            }
        if i.hmm_id == "xJwofaGb0EIZrSxSeZL5xS6thEM7ck7U":
            assert {s: getattr(i, s, None) for s in i.__slots__} == {
                "hmm_id": "xJwofaGb0EIZrSxSeZL5xS6thEM7ck7U",
                "env_from": 95,
                "env_to": 147,
                "seq_pro_score": 45.5,
                "evalue": -8,
                "i_evalue": -8,
                "domain_bias": 0.1,
                "domain_score": 45.5,
                "seq_pro_bias": 0.1,
                "hmm_from": 3,
                "hmm_to": 54,
                "ali_from": 96,
                "ali_to": 147,
                "exponentialized": True,
            }
