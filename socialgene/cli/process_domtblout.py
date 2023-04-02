# python dependencies
import argparse
import csv
from pathlib import Path
from socialgene.base.molbio import Domain

# external dependencies

# internal dependencies
from socialgene.base.socialgene import SocialGene
from socialgene.utils.logging import log
from socialgene.config import env_vars


parser = argparse.ArgumentParser(description="Create *.locus_to_protein a genbank file")
parser.add_argument(
    "--domtblout_file",
    metavar="filepath",
    help="input domtblout file. **Must** be an HMMSEARCH (**not** hmmscan) domtblout file ",
    required=True,
)
parser.add_argument(
    "--outpath",
    metavar="filepath",
    help="Output filepath",
    required=True,
)
parser.add_argument(
    "--sort",
    metavar="boolean",
    help="Read all and sort (to make file reproducible for testing)",
    required=False,
    default=False,
    action=argparse.BooleanOptionalAction,
)


def read_domtblout_write_tsv(domtblout_file, outpath, should_sort):
    if Path(outpath).exists():
        log.error(f"FileExistsError: {outpath}")
        raise FileExistsError
    if not Path(domtblout_file).exists():
        log.error("FileNotFoundError")
        raise FileNotFoundError
    socialgene_object = SocialGene()
    _domain_counter = 0
    if should_sort:
        big_list = list()
        for i in socialgene_object._parse_domtblout(
            input_path=domtblout_file, hmmsearch_or_hmmscan="hmmsearch"
        ):
            domain_obj = Domain(**i)
            _temp = [i["protein_id"]]
            _temp.extend(list(domain_obj.get_dict().values()))
            big_list.add(_temp)
        big_list.sort()
        with open(outpath, "a") as f:
            tsv_writer = csv.writer(f, delimiter="\t")
            for i in socialgene_object._parse_domtblout(
                input_path=domtblout_file, hmmsearch_or_hmmscan="hmmsearch"
            ):
                tsv_writer.writerow(_temp)

    else:
        with open(outpath, "a") as f:
            tsv_writer = csv.writer(f, delimiter="\t")
            socialgene_object = SocialGene()
            for i in socialgene_object._parse_domtblout(
                input_path=domtblout_file, hmmsearch_or_hmmscan="hmmsearch"
            ):
                domain_obj = Domain(**i)
                _temp = [i["protein_id"]]
                _temp.extend(list(domain_obj.get_dict().values()))
                tsv_writer.writerow(_temp)
                _domain_counter += 1
    log.info(f"Wrote {str(_domain_counter)} domains to {outpath}")


def main():
    log.info(f"Socialgene variables: \n{env_vars}")
    # print(Panel(f"Socialgene variables: {env_vars}"))
    args = parser.parse_args()
    log.info(f"sg_process_domtblout will append data to {args.outpath}")
    read_domtblout_write_tsv(domtblout_file=args.domtblout_file, outpath=args.outpath)


if __name__ == "__main__":
    main()
