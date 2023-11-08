import csv
from pathlib import Path

import socialgene.utils.file_handling as fh
from socialgene.parsers.hmmmodel import HmmModel


class HmmInfo:
    def __init__(self, all_hmms_path):
        """
        The function initializes an object with a list of columns and a path to a file.

        Args:
          all_hmms_path: The `all_hmms_path` parameter is a string that represents the path to a file
        containing data for all HMMs.
        """
        self.columns = list(HmmModel().all_attributes().keys())
        self.all_hmms_path = all_hmms_path
        self.all_hmms_data = list()

    def read_tsv(self):
        """
        The function reads a TSV file, validates the number of columns, and appends the data to a list.
        """
        with fh.open_read(self.all_hmms_path) as h:
            for line in h:
                line_vals = [
                    None if v.strip() == '""' else v.strip() for v in line.split("\t")
                ]
                if not len(line_vals) == 16:
                    raise ValueError(f"Expected 16 columns, not {len(line_vals)}")
                if not line_vals == self.columns:
                    # Check for header, if exists, skip reading it in
                    self.all_hmms_data.append(tuple(line_vals))

    def write_nr_hmm_nodes(self, outpath="sg_hmm_nodes"):
        """This writes the HMM node file

        Args:
            single_line_list (list): list representing one row of the HMM metadata file
        """
        with open(outpath, "w") as handle:
            for i in list({i[9] for i in self.all_hmms_data}):
                handle.write(i)
                handle.write("\n")

    def write_hmm_source_nodes(self, outdir="."):
        """Split the HMM metadata into parts based on the database it comes from

        Args:
            single_line_list (list): list representing one row of the HMM metadata file
        """

        for hmm_database_name in {i[0] for i in self.all_hmms_data}:
            # "_hmm_source" is to help glob these files within nextflow
            outpath = Path(outdir, f"{hmm_database_name}_hmm_source")
            with open(outpath, "w") as handle:
                tsv_writer = csv.writer(
                    handle, delimiter="\t", quotechar='"', quoting=csv.QUOTE_MINIMAL
                )
                for i in (i for i in self.all_hmms_data if i[0] == hmm_database_name):
                    tsv_writer.writerow(i[1:])


# write a context manager class that takes a string and returns the string
class PassThrough:
    def __init__(self, string):
        self.string = string

    def __enter__(self):
        return self.string

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass
