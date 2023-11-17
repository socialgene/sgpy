import gzip
import os
import tempfile
from pathlib import Path

from socialgene.cli.nextflow.process_domtblout import main
from socialgene.config import env_vars

FIXTURE_DIR = os.path.dirname(os.path.realpath(__file__))
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)
tbl_path = Path(
    FIXTURE_DIR,
    "data",
    "domtblout",
    "raw",
    "1f97a79b985d9b45262d65113897faee.domtblout.gz",
)
expected_filtered_path = Path(
    FIXTURE_DIR, "data", "domtblout", "processed", "filtered.domtblout.gz"
)
expected_not_filtered_path = Path(
    FIXTURE_DIR, "data", "domtblout", "processed", "not_filtered.domtblout.gz"
)


def test_filtered():
    env_vars["HMMSEARCH_IEVALUE"] = 0.1
    with tempfile.NamedTemporaryFile() as temp_path:
        main(
            args=[
                f"--input={tbl_path}",
                f"--outpath={temp_path.name}",
                "--ievaluefilter",
            ],
        )
        with gzip.open(temp_path.name, "rt") as h:
            test_data = h.readlines()
        with gzip.open(expected_filtered_path, "rt") as h:
            expected_data = h.readlines()
        assert test_data == expected_data


def test_not_filtered():
    env_vars["HMMSEARCH_IEVALUE"] = 0.1
    with tempfile.NamedTemporaryFile() as temp_path:
        main(
            args=[
                f"--input={tbl_path}",
                f"--outpath={temp_path.name}",
            ],
        )
        with gzip.open(temp_path.name, "rt") as h:
            test_data = h.readlines()
        with gzip.open(expected_not_filtered_path, "rt") as h:
            expected_data = h.readlines()
        assert test_data == expected_data


# "--ievaluefilter"
