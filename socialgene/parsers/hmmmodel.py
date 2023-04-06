# python dependencies
from collections import OrderedDict
import csv
from pathlib import Path
import os
import re
from typing import List
from dataclasses import dataclass, field

# external dependencies

# internal dependencies
import socialgene.hashing.hashing as hasher
import socialgene.utils.file_handling as fh
from socialgene.neo4j.schema.define_hmmlist import HMM_SOURCES
from abc import ABC, abstractmethod
from rich import inspect


@dataclass
class HmmModel:
    rel_path: str = None
    model_source: str = None
    hash: str = None
    HMMER3_f: str = None
    NAME: str = None
    ACC: str = None
    DESC: str = None
    LENG: str = None
    MAXL: str = None
    ALPH: str = None
    RF: str = None
    MM: str = None
    CONS: str = None
    CS: str = None
    MAP: str = None
    DATE: str = None
    COM: str = None
    NSEQ: str = None
    EFFN: str = None
    CKSUM: str = None
    GA: str = None
    TC: str = None
    NC: str = None
    STATS: List[str] = field(default_factory=lambda: [])
    HMM: str = None
    COMPO: str = None
    UNKNOWN: List[str] = field(default_factory=lambda: [])
    MODEL: List[str] = field(default_factory=lambda: [])

    def add_attr(self, line):
        """This checks if a line is an HMM attribute (lines above model) and adds them if they are
        Args:
            line (str): string of text coming from file parse/read
        """
        line_val = line.strip().split(maxsplit=1)
        if line.startswith("HMMER3"):
            self.HMMER3_f = line
        elif line_val[0] == "STATS":
            self.STATS.append(line_val[1])
        elif line_val[0] in self.__dataclass_fields__.keys():
            setattr(self, line_val[0], line_val[1])
        else:
            self.UNKNOWN.append(line)

    def add_model(self, line):
        self.MODEL.append(line)

    def add_model_hash(self):
        self.hash = hasher.sha512_hash("".join(self.MODEL))

    def _write_gen_attr_str(self, attr, h):
        if getattr(self, attr):
            n_spaces_offset = len("STATS") + 1 - len(attr)
            h.write(f"{attr}{' ' * n_spaces_offset}{getattr(self, attr)}")
            h.write("\n")

    def write(self, outpath, mode):
        with open(outpath, mode=mode) as h:
            h.write(getattr(self, "HMMER3_f"))
            # write hash if...
            self._write_gen_attr_str("NAME")
            self._write_gen_attr_str("ACC")
            for i in [
                "DESC",
                "LENG",
                "MAXL",
                "ALPH",
                "RF",
                "MM",
                "CONS",
                "CS",
                "MAP",
                "DATE",
                "COM",
                "NSEQ",
                "EFFN",
                "CKSUM",
                "GA",
                "TC",
                "NC",
            ]:
                self._write_gen_attr_str(i)
            for i in self.STATS:
                h.write(f"STATS {i}")
                h.write("\n")
            for i in self.MODEL:
                h.write(i)
            h.write("//")
            h.write("\n")


class HMM:
    def __init__(self) -> None:
        self.models = []
        self.temp_model = HmmModel()

    def read(self, filepath):
        # _temp is used to switch from attribute addition to HMM model additon by switching the function
        # to HMM model's functon after seeing "HMM ", to the end of the model
        # then switch back before the next model starts
        _temp = HmmModel.add_attr
        with fh.open_file(filepath) as h:
            for line in h:
                if line.startswith("HMM "):
                    _temp = HmmModel.add_model
                if line.startswith("//"):
                    self.models.append(self.temp_model)
                    self.temp_model = HmmModel()
                    _temp = HmmModel.add_attr
                _temp(line)

    def write(self, outpath):
        # if only self.read(), then this should write exactly the same file back out
        # (input/output files will have the same hash)
        for i in self.models:
            i.write(outpath, mode="a")
