import argparse
import csv
from multiprocessing import Pool
from pathlib import Path

from socialgene.base.molbio import Domain
from socialgene.base.socialgene import SocialGene
from socialgene.config import env_vars
from socialgene.utils.logging import log

parser = argparse.ArgumentParser(description="Create *.locus_to_protein a genbank file")

parser.add_argument(
    "--input",
    help="One input domtblout file. If a directory, must also use --glob. **Must** be an HMMSEARCH (**not** hmmscan) domtblout file",
    required=True,
    type=str,
    nargs="+",
)
parser.add_argument(
    "--glob",
    help="Glob to input domtblout files. **Must** be an HMMSEARCH (**not** hmmscan) domtblout file ",
    required=False,
    type=str,
)

parser.add_argument(
    "--outpath",
    help="Output filepath",
    required=True,
)
parser.add_argument(
    "--cpus", help="Number of cpus", required=False, default=1, type=int
)


def process_domtblout_file(domtblout_file):
    domain_counter = 0
    socialgene_object = SocialGene()
    _to_return = []
    for i in socialgene_object._parse_domtblout(
        input_path=domtblout_file, hmmsearch_or_hmmscan="hmmsearch"
    ):
        if float(i["i_evalue"]) <= env_vars["HMMSEARCH_IEVALUE"]:
            domain_obj = Domain(**i, exponentialized=True)
            _temp = [i["external_id"]]
            _temp.extend(list(domain_obj.tsv_attributes().values()))
            domain_counter += 1
            _to_return.append(_temp)
    return (_to_return, domain_counter)


def main():
    args = parser.parse_args()
    hmm_files = []
    for i in (Path(i) for i in list(args.input)):
        if i.is_dir():
            if args.glob:
                for ii in list(Path(i).glob(args.glob)):
                    hmm_files.append(ii)
        elif i.is_file():
            hmm_files.append(i)
    if not hmm_files:
        if args.glob:
            log.warning(
                f"No domtblout files found in {args.input} with glob {args.glob}"
            )
        else:
            log.warning(f"No domtblout files found in {args.input}")
    else:
        domain_counter = 0
        log.info(f"Socialgene variables: \n{env_vars}")
        # print(Panel(f"Socialgene variables: {env_vars}"))
        args = parser.parse_args()
        log.info(f"sg_process_domtblout will append data to {args.outpath}")
        p = Pool(args.cpus)
        with open(args.outpath, "w") as f:
            tsv_writer = csv.writer(f, delimiter="\t")
            for result in p.imap(process_domtblout_file, hmm_files):
                for i in result[0]:
                    tsv_writer.writerow(i)
                domain_counter += result[1]
        log.info(f"Wrote {str(domain_counter)} domains to {args.outpath}")


if __name__ == "__main__":
    main()
