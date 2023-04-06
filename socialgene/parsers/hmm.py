# python dependencies
from collections import OrderedDict
import csv
from pathlib import Path
import os
import re

# external dependencies

# internal dependencies
import socialgene.hashing.hashing as hasher
import socialgene.utils.file_handling as fh
from socialgene.neo4j.schema.define_hmmlist import HMM_SOURCES
from abc import ABC, abstractmethod


class HMMCategories:
    """Handles parsing rules specific to certain source databases"""

    def __init__(self):
        self.category = None
        self.subcategory = None

    def parse(self):
        """Extract/Define categories based on parent folder. This is the default scenario"""
        self.category = str(Path(self.rel_path).name)
        self.subcategory = None

    def parse_antismash(self):
        """Extract/Define categories of antismash models (based on antismash directory structure)"""
        # init to None for case of no subcategory
        category = None
        subcategory = None
        temp = str(self.rel_path)
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
                subcategory = temp_split[5].rstrip(".hmm")
        else:
            category = temp_split[2]
        self.category = category
        self.subcategory = subcategory


# The reason we can't just parse everything, hash it and relate back to the names is because
# there is no unifying naming convention/accessions for HMMS across databases. The same "NAME"
# can be located in multiple directories and across databases.
class HMMParser(HMMCategories):
    # the order of hmm_dbs is important because it dictates
    # how the nr model list is created. Later entries override earlier entries
    hmm_dbs = HMM_SOURCES

    def __init__(self):
        super().__init__()
        self.rel_path = None
        self.model_source = None
        self.source_counter = 0
        self.total_hmms_counter = 0
        self.model_text_dict = {}
        self.model_info_dict = {}
        self.single_model_dict = self.blank_model_dict()
        self.temp_list = []
        self.temp_list2 = []
        self.nr_models = {}

    @staticmethod
    def blank_model_dict():
        """Template dictionary to hold the desired info of a single hmm model"""
        return OrderedDict(
            {
                "source": None,
                "rel_path": None,
                "name": None,
                "acc": None,
                "description": None,
                "date": None,
                "sha512t24u": None,
                "model_length": None,
                "category": None,
                "subcategory": None,
            }
        )

    @staticmethod
    def check_if_hmmer3(input_str):
        """Check model to ensure that it has been upconverted to HMMER3

        Args:
            input_str (str): string containing model's hmmer version

        Raises:
            Exception: [description]
        """
        if not input_str.startswith("HMMER3/f"):
            raise Exception("HMMER3 expected")

    def read_body(self, line):
        """Parse the main body of an HMM model

        Args:
            line (str): string of text coming from file parse/read
        """
        if line.startswith("ACC "):
            self.single_model_dict["acc"] = line.strip().split(maxsplit=1)[1]
            # don't include in save model (ie not supposed to have a `temp_list.append(line)`)
        elif line.startswith("NAME "):
            self.single_model_dict["name"] = line.strip().split(maxsplit=1)[1]
        elif line.startswith("DESC "):
            self.single_model_dict["description"] = line.strip().split(maxsplit=1)[1]
            # don't include in save model (ie not supposed to have a `temp_list.append(line)`)
        elif line.startswith("DATE "):
            self.single_model_dict["date"] = line.strip().split(maxsplit=1)[1]
            # don't include in save model (ie not supposed to have a `temp_list.append(line)`)
        else:
            # Create a list of lines that correspond to just the model part of the profile
            # This starts with the line that begins with "HMM " and runs until the line containning "//"
            # this is only used to create a hash to find identical models
            if line.startswith("HMM ") or len(self.temp_list2) > 0:
                self.temp_list2.append(line)
            self.temp_list.append(line)

    def end_of_model(self, line):
        """Actions to perform when reaching the end of parsing an hmm model

        Args:
            line (str): string of text coming from file parse/read
        """
        self.temp_list.append(line)
        self.model_text_dict[self.total_hmms_counter] = self.temp_list
        # temp_list2 contains only the actual hmm model
        self.check_if_hmmer3(self.model_text_dict[self.total_hmms_counter][0])
        self.single_model_dict["sha512t24u"] = hasher.sha512t24u_hasher(
            "".join(self.temp_list2)
        )
        self.single_model_dict["model_length"] = self.temp_list2[-3].split()[0]
        self.single_model_dict["rel_path"] = str(self.rel_path)[
            :-11
        ]  # [:-11] is to remove '_socialgene' suffix
        self.single_model_dict["source"] = str(self.model_source)
        # This changes the model name/id in case of duplicate names within a source database
        self.single_model_dict["name"] = "".join(
            [str(self.single_model_dict["name"]), "_", str(self.source_counter)]
        )
        if self.single_model_dict["source"] == "antismash":
            self.parse_antismash()
        if self.single_model_dict["source"] == "prism":
            self.parse_prism()
        if self.single_model_dict["source"] == "bigslice":
            self.parse_bigslice()
        self.single_model_dict["category"] = self.category
        self.single_model_dict["subcategory"] = self.subcategory
        # Save the parsed hmm info
        self.model_info_dict[self.total_hmms_counter] = self.single_model_dict
        # reset the for-loop dictionary and temp lists
        self.single_model_dict = self.blank_model_dict()
        self.temp_list = []
        self.temp_list2 = []
        self.category = None
        self.subcategory = None
        self.source_counter += 1
        self.total_hmms_counter += 1

    def parse_single_model_file(
        self, input_dir, input_path, model_source, input_rel_path=None
    ):
        """Parse a single file of hmm model(s)

        Args:
            input_dir (Path): Path to top hmm directory, used to help categorize hmm models
            input_path (Path): Path to a single file of hmm models
            model_source (str): source of the models (e.g. 'antismash', 'pfam')
            input_rel_path (str): ignore, for testing
        """
        self.source_counter = 0
        self.model_source = model_source
        input_dir = Path(input_dir)
        input_path = Path(input_path)
        if input_rel_path is not None:
            # used for testing without downloading everything or setting up mock dirs
            self.rel_path = input_rel_path
        else:
            self.rel_path = input_path.relative_to(input_dir)
        with fh.open_file(input_path) as f:
            for line in f:
                if line.strip() == "//":
                    self.end_of_model(line)
                else:
                    self.read_body(line)

    def change_model_name_to_hash(self):
        # For each model, change the 'NAME' line to the model's hash
        for key, model in self.model_text_dict.items():
            # Because the model and parsed model dicts have the same dict keys
            # we pull the model hash value from the same-keyed "model_info_dict" dictionary
            # and replace the NAME with the hash in the model text
            model.insert(
                1,
                "".join(
                    ["NAME  ", self.model_info_dict.get(key).get("sha512t24u"), "\n"]
                ),
            )

    def link_to_latest_pfam(self):
        """Finds latest pfam models by version number and links old models to the new ones
        (e.g. prevent PF02826.20 and PF02826.22 from being included twice)
        """
        # Match 'PF#####'
        re_pfam_broad = re.compile("^PF[0-9]{5,5}")

        # Get all pfam accessions
        pfam_acs1 = {
            v["acc"]
            for v in self.model_info_dict.values()
            if v["acc"] and re_pfam_broad.match(v["acc"])
        }
        # Split versions into list of tuples PFXXXXX.XX -> (PFXXXXX, XX)
        pfam_accs = [(i.split(".")[0], int(i.split(".")[1])) for i in pfam_acs1]
        # Crete dictionary {PFXXXXX: highest_version}
        pfam_version_dict = dict()
        for key, value in pfam_accs:
            if key in pfam_version_dict:
                pfam_version_dict[key] = max(value, pfam_version_dict[key])
            else:
                pfam_version_dict[key] = value

        # find model hashes for highest pfam version
        # Loop through HMM models and replace versioned PFAMs' hash with the highest version's hash
        # add message that it's linked to most recent to acc
        highest_pfam_accs = [f"{k}.{v}" for k, v in pfam_version_dict.items()]
        highest_pfam_accs_shas = {
            v["acc"]: v["sha512t24u"]
            for v in self.model_info_dict.values()
            if v["acc"] in highest_pfam_accs
        }
        for k, v in self.model_info_dict.items():
            if v["acc"] and re_pfam_broad.match(v["acc"]):
                temp_split = v["acc"].split(".")
                if int(temp_split[1]) < pfam_version_dict[temp_split[0]]:
                    pf_id = temp_split[0]
                    highest_version = pfam_version_dict[pf_id]
                    v["acc"] = f"was_{v['acc']}_used_{pf_id}.{highest_version}"
                    v["sha512t24u"] = highest_pfam_accs_shas[
                        f"{pf_id}.{highest_version}"
                    ]

    def create_nr_hmm_dict(self):
        self.nr_models = {}
        self.link_to_latest_pfam()
        for k, v in self.model_info_dict.items():
            if self.model_info_dict[k].get("source") in self.hmm_dbs:
                if isinstance(
                    self.model_info_dict[k].get("acc"), str
                ) and self.model_info_dict[k]["acc"].startswith("was_"):
                    # Skip if it is an outdated pfam
                    continue
                else:
                    self.nr_models[
                        self.model_info_dict[k].get("sha512t24u")
                    ] = self.model_text_dict[k]

    def write_hmms(self, outdir: str = None, n_files: int = 1):
        """Write a non-redundant HMM file, optionally split into n-files"""
        n_models = len(self.nr_models)
        n_files = int(n_files)
        if n_models <= n_files:
            n_per_iter = n_files
        else:
            n_per_iter = int(n_models / n_files)
        written_hmm_counter = 0
        file_counter = 1
        iter_counter = 0
        hmm_filename = f"socialgene_nr_hmms_file_{file_counter}_of_{n_files}.hmm"
        while written_hmm_counter < n_models:
            for hmm in self.nr_models.values():
                if iter_counter >= n_per_iter and file_counter < n_files:
                    file_counter += 1
                    iter_counter = 0
                    hmm_filename = (
                        f"socialgene_nr_hmms_file_{file_counter}_of_{n_files}.hmm"
                    )
                # TODO: bring out of for loop
                with open(os.path.join(outdir, hmm_filename), "a") as hmm_file:
                    for hmm_line in hmm:
                        _ = hmm_file.write(hmm_line)
                iter_counter += 1
                written_hmm_counter += 1

    def write_all_hmm_info(self, outdir):
        with open(os.path.join(outdir, "all_hmms.tsv"), "w") as tsv_file:
            all_hmms_file_writer = csv.DictWriter(
                tsv_file, self.blank_model_dict(), delimiter="\t"
            )
            all_hmms_file_writer.writeheader()
            for line in list(self.model_info_dict.values()):
                _ = all_hmms_file_writer.writerow(line)
