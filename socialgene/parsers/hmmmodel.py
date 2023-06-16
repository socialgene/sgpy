from typing import List

import csv
import os
import re
from collections import OrderedDict, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import socialgene.utils.file_handling as fh
from socialgene.hashing.hashing import hasher
from socialgene.utils.logging import log

re_pfam_broad = re.compile("^PF[0-9]{5,5}")


HMM_SOURCES = [
    "amrfinder",
    "antismash",
    "bigslice",
    "classiphage",
    "ipresto",
    "local",
    "pfam",
    "prism",
    "resfams",
    "tigrfam",
    "virus_orthologous_groups",
]


def parse_hmmlist_input(input):
    # Filter hmm databases based on input list of hmm database names or "all"
    # accept "all" as a list or string
    if isinstance(input, str):
        input = [input]
    if "all" in input:
        _hmms = HMM_SOURCES
    else:
        _hmms = [i for i in input if i in HMM_SOURCES]
    return _hmms


@dataclass
class HmmModel:
    _n: int = None
    _base_dir: str = None
    _abs_path: str = None
    _rel_path: str = None
    _model_source: str = None
    _category: str = None
    _subcategory: str = None
    _hash: str = None
    _new_hash: str = None
    _pfam_accession: str = None
    _pfam_version: int = 0
    _notes: List[str] = field(default_factory=lambda: [])
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
        stripped_line = line.strip()
        line_val = stripped_line.split(maxsplit=1)
        line_val = [i.rstrip() for i in line_val]
        if stripped_line.startswith("HMMER3"):
            self.HMMER3_f = stripped_line
        elif line_val[0] == "STATS":
            self.STATS.append(line_val[1])
        elif line_val[0] in self.__dataclass_fields__.keys():
            setattr(self, line_val[0], line_val[1])
        else:
            self._unknown.append(stripped_line)

    def add_model(self, line):
        self.MODEL.append(line.rstrip())

    def add_model_hash(self):
        self._hash = hasher("".join(self.MODEL))
        self._new_hash = self._hash

    def find_pfam_accessions(self):
        try:
            if re_pfam_broad.match(self.ACC):
                temp = self.ACC.split(".", maxsplit=1)
                self._pfam_accession = temp[0]
                self._pfam_version = int(temp[1])
        except Exception:
            # no error, just don't assign pfam vars
            pass

    def _get_gen_attr_str(self, attr):
        if getattr(self, attr):
            n_spaces_offset = len("STATS") + 1 - len(attr)
            return f"{attr}{' ' * n_spaces_offset}{getattr(self, attr)}"

    def _get_name_as_hash(self):
        attr = "NAME"
        n_spaces_offset = len("STATS") + 1 - len(attr)
        return f"{attr}{' ' * n_spaces_offset}{getattr(self, '_new_hash')}"

    def model_string(self, hash_as_name=False):
        res = []
        # this is spread across multiple calls because likely to include
        # intermediate if/else statemets later
        res.append(getattr(self, "HMMER3_f"))
        # write hash if...
        if hash_as_name:
            res.append(self._get_name_as_hash())
        else:
            res.append(self._get_gen_attr_str("NAME"))
        for i in [
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
            # not all attr are required, check if they exist
            if hasattr(self, i) and getattr(self, i):
                res.append(self._get_gen_attr_str(i))
        for i in self.STATS:
            res.append(f"STATS {i}")
        # write any unknown headers that were present in the original model description
        if getattr(self, "_unknown"):
            for i in getattr(self, "_unknown"):
                res.append(i)
        # write model
        for i in self.MODEL:
            res.append(i)
        res.append("//")
        return "\n".join([i for i in res])

    def write(self, outpath, mode, hash_as_name=False):
        with open(outpath, mode=mode) as h:
            h.write(self.model_string(hash_as_name=hash_as_name))
            h.write("\n")

    def _assign_relative_path(self):
        input_dir = Path(self._base_dir)
        input_path = Path(self._abs_path)
        self._rel_path = input_path.relative_to(input_dir)

    def assign_category(self, input_rel_path=None):
        """Parse antismash categories using antismash directory structure

        Args:
            input_dir (Path): Path to top of hmm directory, used to help categorize hmm models
            input_path (Path): Path to a single file of hmm model(s)
            input_rel_path (str): ignore, for testing
        """
        if not self._base_dir:
            return
        input_dir = Path(self._base_dir)
        if input_rel_path is not None:
            # used for testing without downloading everything or setting up mock dirs
            # TODO: mock instead
            rel_path = input_rel_path
        else:
            self._assign_relative_path()
            rel_path = self._rel_path
        self._model_source = Path(input_dir).name
        if self._model_source == "antismash":
            self._antismash_categories(rel_path)
        else:
            self._category
            self._subcategory

    def _antismash_categories(self, rel_path):
        """Extract/Define categories of antismash models (based on antismash directory structure)"""
        category = None
        subcategory = None
        temp = str(rel_path)
        temp_split = temp.split("/")
        if temp_split[2] == "modules":
            if temp_split[3] == "nrps_pks":
                category = temp_split[3]
                if len(temp_split) > 6:
                    subcategory = temp_split[6]
            elif temp_split[3] in ["sactipeptides", "thiopeptides", "lanthipeptides"]:
                if len(temp_split) == 7:
                    category = temp_split[3]
                    subcategory = temp_split[5]
                else:
                    category = temp_split[3]
            else:
                category = temp_split[3]
        elif temp_split[2] == "detection":
            if temp_split[3] == "hmm_detection":
                category = "hmm_detection"
            elif temp_split[3] == "genefunctions":
                category = "genefunctions"
            elif temp_split[3] == "nrps_pks_domains":
                category = "nrps_pks_domains"
                subcategory = temp_split[5].split(".", maxsplit=1)[0]
        else:
            category = temp_split[2]
        self._category = category
        self._subcategory = subcategory

    def _tsv_dict(self):
        """Template dictionary to hold the desired info of a single hmm model"""
        if isinstance(self._notes, list):
            notes = "; ".join(self._notes)
        else:
            notes = str(self._notes)

        if self._new_hash:
            hash_to_use = self._new_hash
        else:
            hash_to_use = self._hash

        return OrderedDict(
            {
                "id": f"{self._model_source}_{self._n}",
                "source": self._model_source,
                "rel_path": self._rel_path,
                "name": self.NAME,
                "acc": self.ACC,
                "notes": notes,
                "description": self.DESC,
                "date": self.DATE,
                "hash": self._hash,
                "hash_used": hash_to_use,
                "model_length": len(self.MODEL),
                "category": self._category,
                "subcategory": self._subcategory,
                "ga": self.GA,
                "tc": self.TC,
                "nc": self.NC,
            }
        )


class HmmParse:
    def __init__(self) -> None:
        self.models = dict()
        self.dict_key_ind = 0
        self.temp_model = HmmModel()

    def read(self, filepath, base_dir=None):
        # _temp is used to switch from attribute addition to HMM model additon by switching the function
        # to HMM model's functon after seeing "HMM ", to the end of the model
        # then switch back before the next model starts
        _temp = self.temp_model.add_attr
        with fh.open_file(filepath) as h:
            self.temp_model._base_dir = base_dir
            self.temp_model._abs_path = filepath
            self.temp_model._assign_relative_path()
            self.temp_model.assign_category()
            for line in h:
                if line.startswith("HMM "):
                    _temp = self.temp_model.add_model
                if line.startswith("//"):
                    self.temp_model.add_model_hash()
                    self.temp_model.find_pfam_accessions()
                    self.models[self.dict_key_ind] = self.temp_model
                    self.temp_model._n = self.dict_key_ind
                    self.dict_key_ind += 1
                    self.temp_model = HmmModel()
                    self.temp_model._base_dir = base_dir
                    self.temp_model._abs_path = filepath
                    self.temp_model._assign_relative_path()
                    self.temp_model.assign_category()
                    _temp = self.temp_model.add_attr
                    continue
                _temp(line)

    def write_all(self, outpath, hash_as_name=False):
        # if only self.read(), then this should write exactly the same file back out
        # (input/output files will have the same hash)
        for i in self.models.values():
            _ = i.write(outpath, mode="w", hash_as_name=hash_as_name)

    def write_metadata_tsv(self, outdir, header=False):
        with open(os.path.join(outdir, "all.hmminfo"), "w") as tsv_file:
            all_hmms_file_writer = csv.DictWriter(
                tsv_file, self.temp_model._tsv_dict().keys(), delimiter="\t"
            )
            if header:
                all_hmms_file_writer.writeheader()
            for model in self.models.values():
                _ = all_hmms_file_writer.writerow(model._tsv_dict())

    def write_hmm_node_tsv(self, outdir):
        with open(os.path.join(outdir, "sg_hmm_nodes"), "w") as h:
            hashgen = list({i._new_hash for i in self.models.values()})
            hashgen.sort()
            for i in hashgen:
                h.write(i)
                h.write("\n")


class HmmModelHandler(HmmParse):
    def __init__(self) -> None:
        super().__init__()
        self.cull_index = []
        self._pfam_removed = {}
        self._same_hash_removed = {}

    def get_all_model_hashes(self):
        return [getattr(i, "_hash") for i in self.models.values()]

    def hydrate_cull(self):
        self.cull_index = set(self.models.keys())

    def remove_duplicate_and_old_pfam(self):
        cull_at_start = len(self.cull_index)
        # create a dict with
        redundant_pfam = defaultdict(list)
        for k, v in self.models.items():
            if v._pfam_accession:
                redundant_pfam[v._pfam_accession].append(k)
        redundant_pfam = {k: v for k, v in redundant_pfam.items() if len(v) > 1}
        for pfam_id, model_index_list in redundant_pfam.items():
            z = {i: self.models.get(i)._pfam_version for i in model_index_list}
            max_ind = max(z, key=z.get)
            # remove the key of the max
            max_ind_version = z.pop(max_ind)
            deprecated_ind = z.keys()
            for model_index in deprecated_ind:
                # remove from culled list
                try:
                    self.cull_index.remove(model_index)
                except Exception:
                    pass
                # if the version is not the same as the highest version
                # then change metadata for Neo4j

                if self.models[model_index].ACC != self.models[max_ind].ACC:
                    # replace hash in deprecated models so the model metadata
                    # points to the max_ind model in Neo4j
                    self.models[model_index]._new_hash = self.models[max_ind]._hash
                    # change deprecated model ACC to reflect that an updated model was used
                    self.models[model_index]._notes.append(
                        f"{self.models[model_index].ACC}_updated_to_{self.models[max_ind].ACC}"
                    )
            self._pfam_removed[pfam_id] = {
                "kept": f"{pfam_id}.{max_ind_version}",
                "removed": [f"{pfam_id}.{i}" for i in z.values()],
            }
        log.info(
            f"Updated {cull_at_start - len(self.cull_index)} pfam models to a more recent version"
        )

    def remove_duplicate_hash(self):
        cull_at_start = len(self.cull_index)
        redundant_hash = defaultdict(list)
        for k, v in self.models.items():
            redundant_hash[v._hash].append(k)
        redundant_hash = {k: v for k, v in redundant_hash.items() if len(v) > 1}
        for hash_id, model_index_list in redundant_hash.items():
            pass
            keep_this_one = model_index_list.pop(0)
            for i in model_index_list:
                # might have removed in another process (e.g. pfam deduplication), so check first
                if i in self.cull_index:
                    self.cull_index.remove(i)
                self.models[i]._notes.append(
                    f"{self.models[keep_this_one]._rel_path} was used instead of this model"
                )
            self._same_hash_removed[hash_id] = {
                "removed": model_index_list,
                "used": keep_this_one,
            }
        log.info(
            f"{cull_at_start - len(self.cull_index)} identical HMM models removed from culled list"
        )

    def write(
        self,
        outdir: str,
        n_files: int = 1,
        hash_as_name: bool = False,
        splitcutoffs: bool = False,
    ):
        """Write a non-redundant HMM file, optionally split into n-files, or split based on presence
        of cutoff values in a model

        Args:
            outdir (str): where to write the outfile
            n_files (int, optional): split the files into n-files (if splitcutoffs=true will double this value). Defaults to 1.
            hash_as_name (bool, optional): Replace model NAME field with the hash of model. Defaults to False.
            splitcutoffs (bool, optional): Create two file, one with models that have cutoff values, one with models that don't. Defaults to False.
        """
        if any(Path(outdir).glob("socialgene_nr_hmms*")):
            raise FileExistsError(
                f"Found existing socialgene_nr_hmms* file(s) in {outdir}; cowardly not continuing"
            )
        models_without_ga = [k for k, v in self.models.items() if not v.GA]
        models_without_ga = list(
            set(models_without_ga).intersection(set(self.cull_index))
        )
        models_with_ga = [k for k, v in self.models.items() if v.GA]
        models_with_ga = list(set(models_with_ga).intersection(set(self.cull_index)))
        file_counter = 1

        def hmm_filename_with_ga(file_counter, n_files):
            return (
                f"socialgene_nr_hmms_file_with_cutoffs_{file_counter}_of_{n_files}.hmm"
            )

        def hmm_filename_without_ga(file_counter, n_files):
            return f"socialgene_nr_hmms_file_without_cutoffs_{file_counter}_of_{n_files}.hmm"

        for model_index in models_without_ga:
            self.models[model_index].write(
                outpath=Path(outdir, hmm_filename_without_ga(file_counter, n_files)),
                mode="a",
                hash_as_name=hash_as_name,
            )
        for model_index in models_with_ga:
            self.models[model_index].write(
                outpath=Path(outdir, hmm_filename_with_ga(file_counter, n_files)),
                mode="a",
                hash_as_name=hash_as_name,
            )
