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

    def write(self, outpath, mode):
        with open(outpath, mode=mode) as h:
            h.write(getattr(self, "HMMER3_f"))
            for i in [
                "NAME",
                "ACC",
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
                if getattr(self, i):
                    n_spaces_offset = len("STATS") + 1 - len(i)
                    h.write(f"{i}{' ' * n_spaces_offset}{getattr(self, i)}")
                    h.write("\n")
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
        self.dfmlksdl = HmmModel()

    def read(
        self, filepath="/home/chase/Downloads/aaaa/16S_rRNA_NpmA-NCBIFAM.HMM_socialgene"
    ):
        _temp = self.dfmlksdl.add_attr
        with fh.open_file(filepath) as h:
            for line in h:
                if line.startswith("HMM "):
                    _temp = self.dfmlksdl.add_model
                if line.startswith("//"):
                    self.models.append(self.dfmlksdl)
                    self.dfmlksdl = HmmModel()
                    _temp = self.dfmlksdl.add_attr
                _temp(line)

    def write(self, outpath):
        for i in self.models:
            i.write(outpath, mode="a")
