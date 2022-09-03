#!/usr/bin/env python

# python dependencies
from pathlib import Path
import subprocess

# external dependencies

# internal dependencies
from socialgene.config import env_vars
from socialgene.utils.logging import log
from socialgene.utils.run_subprocess import run_subprocess


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
            log.critical(message)
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
        else:
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
    F1: float = env_vars["HMMSEARCH_F1"],
    F2: float = env_vars["HMMSEARCH_F2"],
    F3: float = env_vars["HMMSEARCH_F3"],
    domE=env_vars["HMMSEARCH_DOME"],
    incE=env_vars["HMMSEARCH_INCE"],
    incdomE=env_vars["HMMSEARCH_INCDOME"],
    E=env_vars["HMMSEARCH_E"],
    Z=env_vars["HMMSEARCH_Z"],
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
        Z,
        "-E",
        E,
        "--incE",
        incE,
        "--incdomE",
        incdomE,
        "--domE",
        domE,
        "--F1",
        F1,
        "--F2",
        F2,
        "--F3",
        F3,
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


# TODO: calls to annotate_with_hmmscan should warn if no hmmpressed files
