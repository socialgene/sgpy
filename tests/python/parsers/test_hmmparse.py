import os
import tempfile
from pathlib import Path

from socialgene.parsers.hmmmodel import HmmModel, HmmParse

tmpdir = tempfile.TemporaryDirectory()


FIXTURE_DIR = os.path.dirname(os.path.realpath(__file__))
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)
FIXTURE_DIR = Path(FIXTURE_DIR, "data", "test_genomes")
gbk_path = Path(FIXTURE_DIR, "lagriamide_mibig_bgc0001946.gbk")


FIXTURE_DIR = os.path.dirname(os.path.realpath(__file__))
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)
hmm_dir = Path(FIXTURE_DIR, "data", "hmms")


def test_read():
    hmm_path = Path(hmm_dir, "pks.hmm")
    a = HmmParse()
    a.read(hmm_path, hmm_path.parent)
    del a.models
    {
        "dict_key_ind": 26,
        "temp_model": HmmModel(
            _n=None,
            _base_dir=hmm_dir,
            _abs_path=hmm_path,
            _rel_path=Path("pks.hmm"),
            _model_source="data",
            _category=None,
            _subcategory=None,
            _hash=None,
            _new_hash=None,
            _pfam_accession=None,
            _pfam_version=0,
            _notes=[],
            _unknown=[],
            HMMER3_f=None,
            NAME=None,
            ACC=None,
            DESC=None,
            LENG=None,
            MAXL=None,
            ALPH=None,
            RF=None,
            MM=None,
            CONS=None,
            CS=None,
            MAP=None,
            DATE=None,
            COM=None,
            NSEQ=None,
            EFFN=None,
            CKSUM=None,
            GA=None,
            TC=None,
            NC=None,
            STATS=[],
            HMM=None,
            COMPO=None,
            MODEL=[],
        ),
    }


def check_read_models():
    hmm_path = Path(hmm_dir, "pks.hmm")
    a = HmmParse()
    a.read(hmm_path, hmm_path.parent)
    # check models are held in a dict and they are all models
    assert list(a.models.keys()) == list(range(0, 26))
    assert isinstance(a.models, dict)
    assert len(a.models) == 26
    assert all([type(i) == HmmModel for i in a.models.values()])


def check_write_models():
    hmm_path = Path(hmm_dir, "single_hmm.hmm")
    a = HmmParse()
    a.read(hmm_path, hmm_path.parent)
    with tempfile.NamedTemporaryFile() as fp:
        a.write_all(outpath=fp.name, hash_as_name=False)
        with open(fp.name, "r") as target:
            with open(hmm_path, "r") as expected:
                assert expected.readlines()[1:] == target.readlines()
