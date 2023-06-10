import argparse
from pathlib import Path

from socialgene.parsers.ncbi_taxonomy import (
    merge_taxonomy,
    process_names_dmp,
    process_nodes_dmp,
)
from socialgene.utils.pandas_utils import write_tsv
from socialgene.utils.untargz import untargz

# tax_url = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
# tax_dump_hash_url = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz.md5"


parser = argparse.ArgumentParser(description="Parse NcbiAssemblies taxonomy")
parser.add_argument(
    "--taxdump_path",
    metavar="filepath",
    help="path to ncbi's taxdump.tar.gz",
    required=True,
)


def main():
    args = parser.parse_args()
    taxdump_path = Path(args.taxdump_path)  # taxdump.tar.gz

    names_handle = untargz(taxdump_path, file_name="names.dmp")
    nodes_handle = untargz(taxdump_path, file_name="nodes.dmp")
    names_df = process_names_dmp(file_handle=names_handle)
    nodes_df = process_nodes_dmp(file_handle=nodes_handle)

    write_tsv(
        input_data=merge_taxonomy(names_df=names_df, nodes_df=nodes_df),
        filepath="nodes_taxid",
    )

    write_tsv(input_data=nodes_df.iloc[:, 0:2], filepath="taxid_to_taxid")


if __name__ == "__main__":
    main()
