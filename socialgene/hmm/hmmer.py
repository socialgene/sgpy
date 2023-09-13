import re
import subprocess
from pathlib import Path
from shutil import which
from typing import Generator

import socialgene.utils.file_handling as fh
from socialgene.config import env_vars
from socialgene.utils.logging import log
from socialgene.utils.run_subprocess import run_subprocess


class HMMER:
    """
    The `HMMER` class provides methods for running HMMER programs, such as `hmmpress` and `hmmscan`, and
    checking for the presence of HMMER installation.
    """

    HMMPRESS_OUTPUT_SUFFIXES = ("h3f", "h3i", "h3m", "h3p")

    def __init__(self):
        self._check_hmmer()
        self.hmm_filepaths_with_cutoffs = []
        self.hmm_filepaths_without_cutoffs = []
        self.hmm_filepaths_other = []

    @staticmethod
    def _check_hmmer():
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

    def _append_hmmpress_suffixes(
        self,
        hmm_filepath,
    ) -> Generator[str, None, None]:
        """Appends each of the hmmpress file extensions to the given filepath
        Returns:
            Generator: new paths of expected hmmpress output files
        """
        return (f"{hmm_filepath.suffix}.{i}" for i in self.HMMPRESS_OUTPUT_SUFFIXES)

    def _check_hmmpress_files_exist(self, hmm_filepath):
        """Checks if all the expected files are present after running hmmpress
        Raises:
            FileNotFoundError: Missing hmmpress input/output files
        """
        expected_files = [
            Path(hmm_filepath.with_suffix(i))
            for i in self._append_hmmpress_suffixes(hmm_filepath)
        ]
        if not all([i.exists() for i in expected_files]):
            return False
        else:
            return True

    def hmmpress(self, hmm_path, force=False):
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
        hmm_path = Path(hmm_path)

        if not str(hmm_path).endswith((".hmm", ".hmm.gz")):
            raise ValueError(
                f"HMM file must have an '.hmm' or '.hmm.gz' extension: {hmm_path}"
            )

        if not hmm_path.exists() and hmm_path.suffix == ".gz":
            if hmm_path.with_suffix("").exists():
                hmm_path = hmm_path.with_suffix("")
            else:
                raise FileExistsError

        if fh.is_compressed(hmm_path).name == "gzip":
            # .with_suffix("") removes .gz but leaves .hmm
            self.hmm_filepath = hmm_path.with_suffix("")
        else:
            self.hmm_filepath = hmm_path

        if not force and self._check_hmmpress_files_exist(hmm_path):
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
        print(command_list)
        run_subprocess(
            command_list=command_list,
        )
        if not self._check_hmmpress_files_exist(hmm_path):
            raise FileNotFoundError("Didn't find expected files after running hmmpress")

    def hmmscan(self, hmm_directory, outdirectory, **kwargs):
        for x in Path(hmm_directory).iterdir():
            if re.search(r"\.hmm$|hmm\.gz$", str(x)):
                if "with_cutoffs" in x.stem:
                    self.hmm_filepaths_with_cutoffs.append(x)
                elif "without_cutoffs" in x.stem:
                    self.hmm_filepaths_without_cutoffs.append(x)
                else:
                    self.hmm_filepaths_other.append(x)
        if not any(
            [
                self.hmm_filepaths_with_cutoffs,
                self.hmm_filepaths_without_cutoffs,
                self.hmm_filepaths_other,
            ]
        ):
            raise FileNotFoundError("No HMM model files found")
        for i in self.hmm_filepaths_with_cutoffs:
            self.hmmpress(i)
            temp = dict(**kwargs) | {
                "use_ga_cutoffs": True,
                "domtblout_path": Path(outdirectory, "with_cutoffs.domtblout"),
            }
            self._hmmscan(i, **temp)
        for i in self.hmm_filepaths_without_cutoffs:
            self.hmmpress(i)
            temp = dict(**kwargs) | {
                "use_ga_cutoffs": False,
                "domtblout_path": Path(outdirectory, "without_cutoffs.domtblout"),
            }
            self._hmmscan(i, **temp)
        for i in self.hmm_filepaths_other:
            self.hmmpress(i)
            temp = dict(**kwargs) | {
                "use_ga_cutoffs": False,
                "domtblout_path": Path(outdirectory, "other.domtblout"),
            }
            self._hmmscan(i, **temp)

    @staticmethod
    def _hmmscan(
        hmm_filepath: str,
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
        use_ga_cutoffs: bool = False,
    ) -> None:
        """
        The `hmmscan` function runs HMMER's hmmscan tool to search for protein domains in a given FASTA
        file using a specified HMM profile, and saves the results in a specified output file.

        Args:
          fasta_path (str): The path to the FASTA file that contains the sequences to be searched
        against the HMM profiles.
          input (str): stdin
          hmm_filepath (str): The `hmm_filepath` parameter is the path to the HMM file that will be used
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
        hmm_filepath = Path(hmm_filepath)
        if input:
            fasta_path = "-"
        if not hmm_filepath.exists():
            raise FileNotFoundError(f"No file found at {hmm_filepath}")
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
            hmm_filepath,
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
