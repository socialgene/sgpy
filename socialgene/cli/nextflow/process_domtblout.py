import argparse
import csv
import gzip
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
    "--ievaluefilter",
    help="If true, will filter hmmsearch results by env['HMMSEARCH_IEVALUE'] default limit set in common_parameters.env.",
    required=False,
    default=None,
    action=argparse.BooleanOptionalAction,
)


def process_domtblout_file(domtblout_file, ievaluefilter, domain_counter):
    socialgene_object = SocialGene()
    for i in socialgene_object._parse_domtblout(
        input_path=domtblout_file, hmmsearch_or_hmmscan="hmmsearch"
    ):
        if ievaluefilter:
            filter_value = float(env_vars["HMMSEARCH_IEVALUE"])
        else:
            # if not set, use infinity (keep everything)
            filter_value = float("inf")
        if float(i["i_evalue"]) <= filter_value:
            domain_obj = Domain(**i, exponentialized=True)
            _temp = [i["external_id"]]
            _temp.extend(list(domain_obj.tsv_attributes().values()))
            domain_counter += 1
            yield _temp


def main(args=None):
    args = parser.parse_args(args)
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
        log.info(f"sg_process_domtblout will write data to {args.outpath}")
        with gzip.open(args.outpath, "wt", compresslevel=3) as f:
            tsv_writer = csv.writer(f, delimiter="\t")
            for domtblout_file in hmm_files:
                for i in process_domtblout_file(
                    domtblout_file,
                    ievaluefilter=args.ievaluefilter,
                    domain_counter=domain_counter,
                ):
                    tsv_writer.writerow(i)
                    domain_counter += 1
        log.info(f"Wrote {str(domain_counter)} domains to {args.outpath}")


if __name__ == "__main__":
    main()
