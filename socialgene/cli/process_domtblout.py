# python dependencies
import argparse
from pathlib import Path

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


def read_domtblout_write_tsv(domtblout_file, outpath):
    if Path(outpath).exists():
        log.error(f"FileExistsError: {outpath}")
        raise FileExistsError
    if not Path(domtblout_file).exists():
        log.error("FileNotFoundError")
        raise FileNotFoundError
    socialgene_object = SocialGene()
    socialgene_object.parse_domtblout(
        input_path=domtblout_file, hmmsearch_or_hmmscan="hmmsearch"
    )
    # socialgene_object.filter_domains()
    socialgene_object.export_all_domains_as_tsv(outpath)


def main():
    log.info(f"Socialgene variables: \n{env_vars}")
    # print(Panel(f"Socialgene variables: {env_vars}"))
    args = parser.parse_args()
    read_domtblout_write_tsv(domtblout_file=args.domtblout_file, outpath=args.outpath)


if __name__ == "__main__":
    main()
