import argparse
import csv

from socialgene.base.molbio import Domain
from socialgene.base.socialgene import SocialGene
from socialgene.config import env_vars
from socialgene.utils.logging import log

parser = argparse.ArgumentParser(description="Create *.locus_to_protein a genbank file")
parser.add_argument(
    "--domtblout_file",
    help="input domtblout file. **Must** be an HMMSEARCH (**not** hmmscan) domtblout file ",
    required=True,
)
parser.add_argument(
    "--outpath",
    help="Output filepath",
    required=True,
)


def read_domtblout_write_tsv(domtblout_file, outpath):
    socialgene_object = SocialGene()
    _domain_counter = 0
    with open(outpath, "a") as f:
        tsv_writer = csv.writer(f, delimiter="\t")
        socialgene_object = SocialGene()
        for i in socialgene_object._parse_domtblout(
            input_path=domtblout_file, hmmsearch_or_hmmscan="hmmsearch"
        ):
            if float(i["i_evalue"]) <= env_vars["HMMSEARCH_IEVALUE"]:
                domain_obj = Domain(**i, exponentialized=True)
                _temp = [i["protein_id"]]
                _temp.extend(list(domain_obj.all_attributes.values()))
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
