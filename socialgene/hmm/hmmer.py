import subprocess
from pathlib import Path
from shutil import which
from typing import Generator

import socialgene.utils.file_handling as fh
from socialgene.config import env_vars
from socialgene.parsers.hmmmodel import HmmParse
from socialgene.utils.logging import log
from socialgene.utils.run_subprocess import run_subprocess

HMMPRESS_OUTPUT_SUFFIXES = ("h3f", "h3i", "h3m", "h3p")


class HMMER:
    """
    The `HMMER` class provides methods for running HMMER programs, such as `hmmpress` and `hmmscan`, and
    checking for the presence of HMMER installation.
    """

    def __init__(self, hmm_filepath=None):
        self._check_hmmer_exists()
        self.input_path = Path(hmm_filepath)
        self.decompressed_hmm_path = None

    @staticmethod
    def _check_hmmer_exists():
        """
        The function checks if HMMER is installed and provides the location if found, otherwise it
        displays a warning message.
        """
        if which("hmmsearch"):
            log.info(
                f"Will be using HMMER programs located in: {Path(which('hmmsearch')).parents[0]}"
            )
        else:
            log.warning(
                "Couldn't find HMMER, make sure it is installed and on your PATH"
            )

    @property
    def has_cutoffs(self):
        # Check if the first model in the file has GAthering cutoffs
        path_to_test = None
        if Path(self.input_path).exists():
            path_to_test = self.input_path
        elif Path(self.decompressed_hmm_path).exists():
            path_to_test = self.decompressed_hmm_path
        else:
            raise FileNotFoundError(
                f"Couldn't find hmm file: {self.input_path} or {self.decompressed_hmm_path}"
            )
        temp = HmmParse()
        temp_read = temp.read_model_generator(filepath=path_to_test)
        first_model = temp_read.__next__()
        return first_model.has_cutoffs

    def _decompress_hmm_file_if_needed(
        self,
    ):
        # hmmpress requires a decompressed hmm model file
        # This function checks if the decompressed file exists
        if not str(self.input_path).endswith((".hmm", ".hmm.gz")):
            raise ValueError(
                f"HMM file must have an '.hmm' or '.hmm.gz' extension: {self.input_path}"
            )
        is_compressed = fh.is_compressed(self.input_path).name == "gzip"
        if is_compressed and self.input_path.suffix != ".gz":
            raise ValueError(
                f"Compressed hmm file must have a '.gz' extension: {self.input_path}"
            )
        if is_compressed and self.input_path.with_suffix("").exists():
            log.debug(
                f"Decompressed hmm file already exists: {self.input_path.with_suffix('')}"
            )
            return
        if is_compressed:
            self.decompressed_hmm_path = fh.gunzip(self.input_path)
            log.debug(
                f"Decompressed hmm file: {self.input_path} -> {self.decompressed_hmm_path}"
            )
        else:
            self.decompressed_hmm_path = self.input_path

    def _append_hmmpress_suffixes(self) -> Generator[str, None, None]:
        """Appends each of the hmmpress file extensions to the given filepath
        Returns:
            Generator: new paths of expected hmmpress output files
        """
        return (f"{self.decompressed_hmm_path}.{i}" for i in HMMPRESS_OUTPUT_SUFFIXES)

    def _check_hmmpress_files_exist(self) -> bool:
        """Checks if all the expected files are present after running hmmpress
        Returns:
            bool: True if all expected files are present, False otherwise
        """
        if all((Path(i).exists() for i in self._append_hmmpress_suffixes())):
            return True
        else:
            log.debug("Not all expected hmmpress files found")
            log.debug(
                f"Missing: {(i for i in self._append_hmmpress_suffixes() if not i.exists())}"
            )
            return False

    def hmmpress(self, force=False):
        """
        The function `hmmpress` checks if the HMM file has the correct extension, determines the file
        path, and skips the hmmpress step if the necessary files already exist.

        Args:
          hmm_path: The path to the HMM file that needs to be processed.
          force: The `force` parameter is a boolean flag that determines whether to force the execution
        of the `hmmpress` function or not. If `force` is set to `True`, the function will be executed
        regardless of whether the necessary files already exist. Defaults to False

        Returns:
          nothing (None).
        """
        self._decompress_hmm_file_if_needed()
        # `self._decompress_hmm_file_if_needed()` populates `self.decompressed_hmm_path`
        if not force and self._check_hmmpress_files_exist():
            log.info(
                "hmmpress outputs found and hmmpress(force=False), skipping hmmpress"
            )
            return
        if force:
            command_list = ["hmmpress", "-f", str(self.decompressed_hmm_path)]
        else:
            command_list = ["hmmpress", str(self.decompressed_hmm_path)]
        command_list = [str(i) for i in command_list]
        run_subprocess(
            command_list=command_list,
        )
        if not self._check_hmmpress_files_exist():
            raise FileNotFoundError("Didn't find expected files after running hmmpress")

    def hmmscan(
        self,
        domtblout_path: str,
        fasta_path: str,
        input: str = None,
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
        """
        The `hmmscan` function runs HMMER's hmmscan tool to search for protein domains in a given FASTA
        file using a specified HMM profile, and saves the results in a specified output file.

        Args:
          fasta_path (str): The path to the FASTA file that contains the sequences to be searched
        against the HMM profiles.
          input (str): use to pass fasta as stdin
          hmmpath (str): The `hmmpath` parameter is the path to the HMM file that will be used
        for the hmmscan.
          domtblout_path (str): The `domtblout_path` parameter is the path where the output file from
        HMMER's hmmscan will be saved. This file will contain the domain hits found in the input
        sequences.
          f1 (float): The parameter `f1` is a float value representing the gathering threshold for the
        first domain architecture. It is used in the HMMER's hmmscan command to set the threshold for
        the first domain architecture. The default value is obtained from the `env_vars` dictionary with
        the key "H
          f2 (float): The parameter `f2` in the `hmmscan` method is a float that represents the
        gathering threshold for the second domain architecture. It is used as an optional argument in
        the `hmmscan` command.
          f3 (float): The parameter `f3` is a float value that represents the threshold for the third
        stage of the HMMER's hmmscan algorithm. It is used to control the trade-off between sensitivity
        and speed during the search. A higher value of `f3` will increase the speed of the search
          dom_e: The `dom_e` parameter is the E-value threshold for reporting domain hits. It is used to
        control the significance threshold for reporting individual domains in the output.
          inc_e: The `inc_e` parameter is the inclusion threshold for the E-value. It is used to
        determine if a hit is significant enough to be included in the output. Hits with E-values below
        the inclusion threshold will be reported.
          incdom_e: The `incdom_e` parameter is the inclusion threshold for domain E-values in HMMER's
        hmmscan. It is used to control the sensitivity of the search by setting a threshold for the
        significance of domain hits. Any domain hit with an E-value below the `incdom_e` threshold
          e: The parameter `e` in the `hmmscan` method is the E-value threshold for reporting domains.
        It is a floating-point value that determines the significance threshold for the reported
        matches. A lower E-value indicates a more significant match.
          z: The parameter `z` in the `hmmscan` method is used to set the reporting threshold for the
        gathering score. The gathering score is a measure of the significance of a match between a
        sequence and a profile hidden Markov model (HMM). The `z` parameter specifies the number of
          seed (int): The `seed` parameter in the `hmmscan` method is an optional integer parameter that
        specifies the random number seed for the HMMER's hmmscan algorithm.
          cpus (int): The `cpus` parameter specifies the number of CPUs (or cores) to be used for
        running the hmmscan process. Defaults to 1
          overwrite: The `overwrite` parameter is a boolean flag that determines whether to overwrite an
        existing file at the `domtblout_path` location. If `overwrite` is set to `False` and a file
        already exists at `domtblout_path`, a `FileExistsError` will be raised. Defaults to False
        """
        self.hmmpress()

        if self.has_cutoffs:
            # base on the first model in the file
            use_ga_cutoffs = True
        else:
            use_ga_cutoffs = False

        if input:
            # set hmmscan arg to read from stdin
            fasta_path = "-"

        domtblout_path = Path(domtblout_path)
        if domtblout_path.exists() and not overwrite:
            raise FileExistsError(f"File already exists at {domtblout_path}")
        c1 = [
            "hmmscan",
            "--seed",
            int(seed),
            "--cpu",
            int(cpus),
            "-Z",
            z,
        ]

        c3 = [
            "--F1",
            f1,
            "--F2",
            f2,
            "--F3",
            f3,
            "--domtblout",
            domtblout_path,
            self.decompressed_hmm_path,
            fasta_path,
        ]
        if use_ga_cutoffs:
            c2 = ["--cut_ga"]
        else:
            c2 = [
                "-E",
                e,
                "--incE",
                inc_e,
                "--incdomE",
                incdom_e,
                "--domE",
                dom_e,
            ]
        command_list = c1 + c2 + c3
        command_list = [str(i) for i in command_list]
        run_subprocess(
            command_list=command_list,
            input=input,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            shell=False,
            capture_output=False,
        )
