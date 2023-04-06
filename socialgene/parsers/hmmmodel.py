# python dependencies
from collections import OrderedDict
import csv
from pathlib import Path
import os
import re
from typing import List
from dataclasses import dataclass, field
from collections import defaultdict

# external dependencies

# internal dependencies
import socialgene.hashing.hashing as hasher
import socialgene.utils.file_handling as fh
from socialgene.neo4j.schema.define_hmmlist import HMM_SOURCES
from abc import ABC, abstractmethod
from rich import inspect

re_pfam_broad = re.compile("^PF[0-9]{5,5}")


@dataclass
class HmmModel:
    _rel_path: str = None
    _model_source: str = None
    _hash: str = None
    _pfam_accession: str = None
    _pfam_version: int = 0
    _unknown: List[str] = field(default_factory=lambda: [])
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
            self._unknown.append(line)

    def add_model(self, line):
        self.MODEL.append(line)

    def add_model_hash(self):
        self._hash = hasher.sha512t24u_hasher("".join(self.MODEL))

    def find_pfam_accessions(self):
        try:
            if re_pfam_broad.match(self.ACC):
                temp = self.ACC.split(".", maxsplit=1)
                self._pfam_accession = temp[0]
                self._pfam_version = int(temp[1])
        except:
            # no error, just don't assign pfam vars
            pass

    def _write_gen_attr_str(self, attr, h):
        if getattr(self, attr):
            n_spaces_offset = len("STATS") + 1 - len(attr)
            h.write(f"{attr}{' ' * n_spaces_offset}{getattr(self, attr)}")
            h.write("\n")

    def write(self, outpath, mode):
        with open(outpath, mode=mode) as h:
            # write model headers
            # this is spread across multiple calls because likely to include
            # intermediate if/else statemets later
            h.write(getattr(self, "HMMER3_f"))
            # write hash if...
            self._write_gen_attr_str("NAME", h)
            self._write_gen_attr_str("ACC", h)
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
                self._write_gen_attr_str(i, h)
            for i in self.STATS:
                h.write(f"STATS {i}")
                h.write("\n")
            # write any unknown headers that were present in the original model description
            if getattr(self, "_unknown"):
                for i in getattr(self, "_unknown"):
                    h.write(i)
                    h.write("\n")
            # write model
            for i in self.MODEL:
                h.write(i)
            h.write("//")
            h.write("\n")


class HmmParse:
    def __init__(self) -> None:
        self.models = dict()
        self.temp_model = HmmModel()

    def read(self, filepath):
        # _temp is used to switch from attribute addition to HMM model additon by switching the function
        # to HMM model's functon after seeing "HMM ", to the end of the model
        # then switch back before the next model starts
        dict_key_ind = 0
        _temp = self.temp_model.add_attr
        with fh.open_file(filepath) as h:
            for line in h:
                if line.startswith("HMM "):
                    _temp = self.temp_model.add_model
                if line.startswith("//"):
                    self.temp_model.add_model_hash()
                    self.temp_model.find_pfam_accessions()
                    self.models[dict_key_ind] = self.temp_model
                    dict_key_ind += 1
                    self.temp_model = HmmModel()
                    _temp = self.temp_model.add_attr
                    continue
                _temp(line)

    def write_all(self, outpath):
        # if only self.read(), then this should write exactly the same file back out
        # (input/output files will have the same hash)
        for i in self.models.values():
            i.write(outpath, mode="a")


class H(HmmParse):
    def __init__(self) -> None:
        super().__init__()
        self.cull_index = []
        self._pfam_removed = {}
        self._same_hash_removed = {}

    def get_all_model_hashes(self):
        return [getattr(i, "_hash") for i in self.models.values()]

    def huh(self):
        self.hashlist = self.get_all_model_hashes()

    def hydrate_cull(self):
        self.cull_index = set(self.models.keys())

    def remove_duplicate_and_old_pfam(self):
        redundant_pfam = defaultdict(list)
        for k, v in self.models.items():
            if v._pfam_accession:
                redundant_pfam[v._pfam_accession].append(k)
        redundant_pfam = {k: v for k, v in redundant_pfam.items() if len(v) > 1}
        for pfam_id, model_index_list in redundant_pfam.items():
            z = {i: self.models.get(i)._pfam_version for i in model_index_list}
            max_ind = max(z, key=z.get)
            # remove the key of the max
            _kept = z.pop(max_ind)
            for i in z.keys():
                self.cull_index.remove(i)
            self._pfam_removed[pfam_id] = {
                "kept": f"{pfam_id}.{_kept}",
                "removed": [f"{pfam_id}.{i}" for i in z.values()],
            }

    def remove_duplicate_hash(self):
        redundant_hash = defaultdict(list)
        for k, v in self.models.items():
            redundant_hash[v._hash].append(k)
        redundant_hash = {k: v for k, v in redundant_hash.items() if len(v) > 1}
        for hash_id, model_index_list in redundant_hash.items():
            _kept = model_index_list.pop()
            for i in model_index_list:
                # might have removed in another process (e.g. pfam deduplication), so check first
                if i in self.cull_index:
                    self.cull_index.remove(i)
            self._same_hash_removed[hash_id] = {
                "removed": model_index_list,
            }

    def write_culled(self, outpath):
        # if only self.read(), then this should write exactly the same file back out
        # (input/output files will have the same hash)
        for i in self.cull_index:
            print(i)
            self.models[i].write(outpath, mode="a")
