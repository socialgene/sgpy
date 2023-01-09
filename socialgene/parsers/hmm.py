# python dependencies
import csv
from pathlib import Path
import os

# external dependencies

# internal dependencies
import socialgene.hashing.hashing as hasher
import socialgene.utils.file_handling as fh

# Add new HMM parsers to IndividualHmmDbParsers


class IndividualHmmDbParsers:
    def __init__(self):
        self.category = None
        self.subcategory = None

    def parse_parent_folder_is_name(self):
        """Extract/Define categories based on parent folder"""
        # init to None for case of no subcategory
        subcategory = None
        temp = str(self.rel_path)
        temp_split = temp.split("/")
        category = temp_split[1]
        self.category = category
        self.subcategory = subcategory

    def parse_prism(self):
        self.parse_parent_folder_is_name()

    def parse_bigslice(self):
        self.parse_parent_folder_is_name()

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
class HMMParser(IndividualHmmDbParsers):
    # the order of hmm_dbs is important because it dictates
    # how the nr model list is created. Later entries override earlier entries
    hmm_dbs = [
        "prism",
        "bigslice",
        "antismash",
        "amrfinder",
        "resfams",
        "tigrfam",
        "pfam",
        "classiphage",
        "virus_orthologous_groups",
        "local",
    ]

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
        return {
            "source": None,
            "source_id": None,
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
        self.single_model_dict["sha512t24u"] = hasher.sha512_hash(
            "".join(self.temp_list2)
        )
        self.single_model_dict["model_length"] = self.temp_list2[-3].split()[0]
        self.single_model_dict["rel_path"] = str(self.rel_path)[
            :-11
        ]  # [:-11] is to remove '_socialgene' suffix
        self.single_model_dict["source"] = str(self.model_source)
        self.single_model_dict["source_id"] = "".join(
            [str(self.model_source), "_", str(self.source_counter)]
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

    def create_nr_hmm_dict(self):
        self.nr_models = {}
        for source_db in self.hmm_dbs:
            for k, v in self.model_info_dict.items():
                if self.model_info_dict[k].get("source") == source_db:
                    self.nr_models[
                        self.model_info_dict[k].get("sha512t24u")
                    ] = self.model_text_dict[k]

    def write_hmms(self, outdir: str = None, n_files: int = 1):
        """Write a non-redundant HMM file, optionally split into x files"""
        n_models = len(self.nr_models.values())
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
