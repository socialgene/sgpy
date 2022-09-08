# python dependencies
from pathlib import Path
import subprocess
import gzip
import shutil

# external dependencies

# internal dependencies
from socialgene.config import env_vars
from socialgene.utils.logging import log
from socialgene.utils.run_subprocess import run_subprocess
import socialgene.utils.file_handling as fh


# Check if the hmm file is still compressed (needs to be decompressed for hmmpress)
def check_pressed_files(hmm_filepath, fail_on_missing=False):
    hmm_filepath = Path(hmm_filepath)
    expected_files = [
        Path(hmm_filepath.with_suffix(i))
        for i in [".hmm.h3f", ".hmm.h3i", ".hmm.h3m", ".hmm.h3p"]
    ]
    if not hmm_filepath.exists():
        raise FileNotFoundError(hmm_filepath)
    if not all([i.exists() for i in expected_files]):
        message = f"Missing hmmpress files for: {hmm_filepath}\n HMMER's hmmpress can be run using the hmmpress() function in socialgene.hmm.hmmscan"
        if fail_on_missing:
            log.error(message)
            raise FileNotFoundError(message)
        else:
            log.info(message)


def hmmpress(hmm_filepath, overwrite=False):
    hmm_filepath = Path(hmm_filepath)
    expected_files = [
        Path(hmm_filepath.with_suffix(i))
        for i in [".hmm", ".hmm.h3f", ".hmm.h3i", ".hmm.h3m", ".hmm.h3p"]
    ]
    if overwrite:
        command_list = ["hmmpress", "-f", str(hmm_filepath)]
    else:
        if all([i.exists() for i in expected_files]):
            log.info(
                "hmmpress outputs found and hmmpress(overwrite=False), skipping hmmpress"
            )
            return
        elif any([i.exists() for i in expected_files[1:]]):
            log.info(
                "Partial hmmpress outputs found; which means hmmpress(overwrite=False), will fail.\n Figure out what your deal is, or run again with hmmpress(force=True)"
            )
            raise FileNotFoundError
        #  Decompress first if hmm file is gzipped
        if fh.is_compressed(hmm_filepath).name == "gzip":
            new_hmm_path = Path(hmm_filepath.parents[0], hmm_filepath.stem)
            with gzip.open(hmm_filepath, "rb") as f_in:
                with open(new_hmm_path, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            hmm_filepath = new_hmm_path
        command_list = ["hmmpress", str(hmm_filepath)]
        command_list = [str(i) for i in command_list]
    run_subprocess(
        command_list=command_list,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    check_pressed_files(hmm_filepath=hmm_filepath, fail_on_missing=True)


def run_hmmscan(
    fasta_path: str,
    input: str,
    hmm_filepath: str,
    domtblout_path: str,
    f1: float = env_vars["HMMSEARCH_F1"],
    f2: float = env_vars["HMMSEARCH_F2"],
    f3: float = env_vars["HMMSEARCH_F3"],
    dom_e=env_vars["HMMSEARCH_DOME"],
    inc_e=env_vars["HMMSEARCH_INCE"],
    incdom_e=env_vars["HMMSEARCH_INCDOME"],
    e=env_vars["HMMSEARCH_E"],
    z=env_vars["HMMSEARCH_Z"],
    seed: int = env_vars["HMMSEARCH_SEED"],
    cpus: int = 1,
    overwrite=False,
):
    hmm_filepath = Path(hmm_filepath)
    if not hmm_filepath.exists():
        raise FileNotFoundError(f"No file found at {hmm_filepath}")
    check_pressed_files(hmm_filepath=hmm_filepath, fail_on_missing=True)
    domtblout_path = Path(domtblout_path)
    if domtblout_path.exists() and not overwrite:
        raise FileExistsError(f"File already exists at {domtblout_path}")
    command_list = [
        "hmmscan",
        "--seed",
        int(seed),
        "--cpu",
        int(cpus),
        "-Z",
        z,
        "-E",
        e,
        "--incE",
        inc_e,
        "--incdomE",
        incdom_e,
        "--domE",
        dom_e,
        "--F1",
        f1,
        "--F2",
        f2,
        "--F3",
        f3,
        "--domtblout",
        domtblout_path,
        hmm_filepath,
        fasta_path,
    ]
    command_list = [str(i) for i in command_list]
    run_subprocess(
        command_list=command_list,
        input=input,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        shell=False,
    )
