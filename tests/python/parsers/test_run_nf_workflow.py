import os
import tempfile
from pathlib import Path

import pytest

from socialgene.cli.nextflow.clean_hmms import run_nf_workflow

FIXTURE_DIR = os.path.dirname(os.path.realpath(__file__))
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)
hmm_dirs = Path(FIXTURE_DIR, "data", "hmm_dirs")

FIXTURE_DIR = os.path.dirname(os.path.realpath(__file__))
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)
hmm_dirs_processed = Path(FIXTURE_DIR, "data", "hmm_dirs_processed")


def test_run_nf_workflow_files_exist():
    with tempfile.TemporaryDirectory() as temp_path:
        run_nf_workflow(
            input_dir=hmm_dirs,
            outdir=temp_path,
            input_glob="**/*.socialgene.hmm.gz",
        )
        assert Path(temp_path, "all.hmminfo").exists()
        assert Path(temp_path, "sg_hmm_nodes").exists()
        assert Path(
            temp_path, "socialgene_nr_hmms_file_with_cutoffs_1_of_1.hmm"
        ).exists()
        assert Path(
            temp_path, "socialgene_nr_hmms_file_without_cutoffs_1_of_1.hmm"
        ).exists()


files = [
    "all.hmminfo",
    "sg_hmm_nodes",
    "socialgene_nr_hmms_file_with_cutoffs_1_of_1.hmm",
    "socialgene_nr_hmms_file_without_cutoffs_1_of_1.hmm",
]


@pytest.mark.parametrize("file", files)
def test_run_nf_workflow_sg_hmm_nodes(file):
    with tempfile.TemporaryDirectory() as temp_path:
        run_nf_workflow(
            input_dir=hmm_dirs,
            outdir=temp_path,
            input_glob="**/*.socialgene.hmm.gz",
        )
        with open(Path(temp_path, file)) as f:
            with open(Path(hmm_dirs_processed, file)) as f2:
                assert f.readlines() == f2.readlines()
