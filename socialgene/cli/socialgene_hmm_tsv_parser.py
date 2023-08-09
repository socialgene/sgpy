import argparse
from pathlib import Path

from socialgene.hmm.hmminfo import HmmInfo

parser = argparse.ArgumentParser(description="Merging tsv")
parser.add_argument(
    "--all_hmms_path",
    help=".sg_hmm_nodes filepath",
    required=True,
)


def main():
    args = parser.parse_args()
    all_hmms_path = Path(args.all_hmms_path)
    a = HmmInfo(all_hmms_path)
    a.read_tsv()
    a.write_nr_hmm_nodes(outpath="sg_hmm_nodes")
    a.write_hmm_source_nodes(outdir=".")


if __name__ == "__main__":
    main()
