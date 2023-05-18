# python dependencies
import csv
from pathlib import Path

# external dependencies

# internal dependencies
from socialgene.parsers.hmmmodel import HmmModel
import socialgene.utils.file_handling as fh

HmmModel


class HmmInfo:
    def __init__(self, all_hmms_path):
        self.columns = list(HmmModel()._tsv_dict().keys())
        self.all_hmms_path = all_hmms_path
        self.all_hmms_data = list()

    def read_tsv(self):
        with fh.open_file(self.all_hmms_path) as h:
            for line in h:
                line_vals = [
                    None if v.strip() == '""' else v.strip() for v in line.split("\t")
                ]
                if not len(line_vals) == 10:
                    print(line_vals)
                    raise ValueError(f"Expected 10 columns, not {len(line_vals)}")
                if not line_vals == self.columns:
                    self.all_hmms_data.append(tuple(line_vals))

    def write_nr_hmm_nodes(self, outpath="sg_hmm_nodes"):
        """This writes the HMM node file

        Args:
            single_line_list (list): list representing one row of the HMM metadata file
        """
        with open(outpath, "w") as handle:
            tsv_writer = csv.writer(
                handle, delimiter="\t", quotechar='"', quoting=csv.QUOTE_MINIMAL
            )
            for i in {(i[6], i[7]) for i in self.all_hmms_data}:
                tsv_writer.writerow(i)

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
