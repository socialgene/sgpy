from typing import Generator

import subprocess
from pathlib import Path
from shutil import which

import socialgene.utils.file_handling as fh
from socialgene.config import env_vars
from socialgene.utils.logging import log
from socialgene.utils.run_subprocess import run_subprocess


class HMMER:
    HMMPRESS_OUTPUT_SUFFIXES = ("h3f", "h3i", "h3m", "h3p")

    def __init__(self):
        self._check_hmmer()
        self.hmm_filepath = None

    @staticmethod
    def _check_hmmer():
        if which("hmmsearch"):
            log.info(
                f"Will be using HMMER programs located in: {Path(which('hmmsearch')).parents[0]}"
            )
        else:
            log.warning(
                "Couldn't find HMMER, make sure it is installed and on your PATH"
            )

    def _append_hmmpress_suffixes(
        self,
    ) -> Generator[str, None, None]:
        """Appends each of the hmmpress file extensions to the given filepath
        Returns:
            Generator: new paths of expected hmmpress output files
        """
        return (
            f"{self.hmm_filepath.suffix}.{i}" for i in self.HMMPRESS_OUTPUT_SUFFIXES
        )

    def _check_hmmpress_files_exist(self):
        """Checks if all the expected files are present after running hmmpress
        Raises:
            FileNotFoundError: Missing hmmpress input/output files
        """
        expected_files = [
            Path(self.hmm_filepath.with_suffix(i))
            for i in self._append_hmmpress_suffixes()
        ]
        expected_files = expected_files

        if not all([i.exists() for i in expected_files]):
            return False
        else:
            return True

    def hmmpress(self, hmm_path, force=False):
        hmm_path = Path(hmm_path)

        if not str(hmm_path).endswith((".hmm", ".hmm.gz")):
            raise ValueError("HMM file must have an '.hmm' or '.hmm.gz' extension")

        if fh.is_compressed(hmm_path).name == "gzip":
            # .with_suffix("") removes .gz but leaves .hmm
            self.hmm_filepath = hmm_path.with_suffix("")
        else:
            self.hmm_filepath = hmm_path

        if not force and self._check_hmmpress_files_exist():
            log.info(
                "hmmpress outputs found and hmmpress(force=False), skipping hmmpress"
            )
            return

        if fh.is_compressed(hmm_path).name == "gzip":
            if force or not self.hmm_filepath.exists():
                fh.gunzip(hmm_path)

        if force:
            command_list = ["hmmpress", "-f", str(self.hmm_filepath)]
        else:
            #  Decompress first if hmm file is gzipped
            command_list = ["hmmpress", str(self.hmm_filepath)]
        command_list = [str(i) for i in command_list]
        run_subprocess(
            command_list=command_list,
        )
        if not self._check_hmmpress_files_exist():
            raise FileNotFoundError("Didn't find expected files after running hmmpress")

    @staticmethod
    def hmmscan(
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
    ) -> None:
        """Run HMMER's hmmscan

        Args:
            fasta_path (str): _description_
            input (str): _description_
            hmm_filepath (str): _description_
            domtblout_path (str): _description_
            f1 (float, optional): _description_. Defaults to env_vars["HMMSEARCH_F1"].
            f2 (float, optional): _description_. Defaults to env_vars["HMMSEARCH_F2"].
            f3 (float, optional): _description_. Defaults to env_vars["HMMSEARCH_F3"].
            dom_e (_type_, optional): _description_. Defaults to env_vars["HMMSEARCH_DOME"].
            inc_e (_type_, optional): _description_. Defaults to env_vars["HMMSEARCH_INCE"].
            incdom_e (_type_, optional): _description_. Defaults to env_vars["HMMSEARCH_INCDOME"].
            e (_type_, optional): _description_. Defaults to env_vars["HMMSEARCH_E"].
            z (_type_, optional): _description_. Defaults to env_vars["HMMSEARCH_Z"].
            seed (int, optional): _description_. Defaults to env_vars["HMMSEARCH_SEED"].
            cpus (int, optional): _description_. Defaults to 1.
            overwrite (bool, optional): _description_. Defaults to False.

        Raises:
            FileNotFoundError: _description_
            FileExistsError: _description_
        """
        hmm_filepath = Path(hmm_filepath)
        if not hmm_filepath.exists():
            raise FileNotFoundError(f"No file found at {hmm_filepath}")
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
            capture_output=False,
        )
