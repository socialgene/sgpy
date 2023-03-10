# python dependencies
import csv
from pathlib import Path
from socialgene.parsers.hmm import HMMParser

# external dependencies

# internal dependencies


# all_hmms_path = "/home/chase/Documents/socialgene_data/examplesssss/socialgene_per_run/hmm_cache/all_hmms.tsv"


class HmmInfo:
    def __init__(self, all_hmms_path):
        self.columns = list(HMMParser.blank_model_dict().keys())
        self.all_hmms_path = all_hmms_path
        self.track_hmm_ids = dict()
        self.all_hmms_data = set()

    def read_tsv(self):
        with open(self.all_hmms_path, "r") as h:
            for line in h:
                line_vals = [i.strip() for i in line.split("\t")]
                if not len(line_vals) == 10:
                    print(line_vals)
                    raise ValueError(f"Expected 10 columns, not {len(line_vals)}")
                if not line_vals == self.columns:
                    self.all_hmms_data.add(tuple(line_vals))

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
            outpath = Path(outdir, hmm_database_name)
            with open(outpath, "w") as handle:
                tsv_writer = csv.writer(
                    handle, delimiter="\t", quotechar='"', quoting=csv.QUOTE_MINIMAL
                )
                for i in (i for i in self.all_hmms_data if i[0] == hmm_database_name):
                    tsv_writer.writerow(i[1:])
