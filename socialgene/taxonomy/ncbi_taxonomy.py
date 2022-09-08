# python dependencies
from pathlib import Path
import argparse

# external dependencies
import pandas as pd

# internal dependencies
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


def process_names_dmp(file_handle):
    df = pd.read_csv(
        file_handle,
        sep="|",
        header=None,
        index_col=False,
        names=[
            "tax_id",
            "name_txt",
            "unique_name",
            "name_class",
        ],
        engine="c",
    )
    # remove extra \t
    df["name_txt"] = df["name_txt"].str.strip()
    df["unique_name"] = df["unique_name"].str.strip()
    df["name_class"] = df["name_class"].str.strip()
    # remove synonyms
    df = df[df["name_class"] == "scientific name"]
    df = df.iloc[:, [0, 1]]
    return df


def process_nodes_dmp(file_handle):
    df2 = pd.read_csv(
        file_handle,
        sep="|",
        header=None,
        index_col=False,
        names=[
            "tax_id",
            "parent_tax_id",
            "rank",
            "embl_code",
            "division_id",
            "inherited_div_flag",
            "genetic_code_id",
            "inherited_GC__flag",
            "mitochondrial_genetic_code_id",
            "inherited_MGC_flag",
            "GenBank_hidden_flag",
            "hidden_subtree_root_flag",
            "comments",
        ],
        engine="c",
    )
    df2["rank"] = df2["rank"].str.strip()
    df2["embl_code"] = df2["embl_code"].str.strip()
    df2["comments"] = df2["comments"].str.strip()
    return df2


def merge_taxonomy(names_df, nodes_df):
    return names_df.merge(nodes_df.iloc[:, [0, 2]], how="inner", on="tax_id")


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
