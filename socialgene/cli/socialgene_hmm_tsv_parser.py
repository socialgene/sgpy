# python dependencies
import argparse
from pathlib import Path

# external dependencies
import pandas as pd

# internal dependencies
from socialgene.utils.pandas_utils import write_tsv
from socialgene.utils.logging import log
from socialgene.neo4j.schema.define_hmmlist import HMM_SOURCES

parser = argparse.ArgumentParser(description="Merging tsv")
parser.add_argument(
    "--all_hmms_path",
    metavar="filepath",
    help="all_hmms.tsv filepath",
    required=True,
)

# socialgene creates three tsv files when creating unique hmm ids
# - 'all_hmms.tsv' contains all the hmm files that it found
#     - some of these are duplicates and so the first column `sha512t24u` is *not* unique
# - 'unique_hmms.tsv' are the hmm models that were passed into the socialgene pipeline to annotate fasta files
#     - the first column `sha512t24u` *is* unique
# - 'removed_hmms.tsv' contains the hmms removed from 'all_hmms.tsv' to create the set of 'unique_hmms.tsv'


def unique_hmm_names(dict_of_dfs):
    for single_df in dict_of_dfs.values():
        input_column = single_df["name"]
        temp_dict = dict(zip(input_column, [0] * len(input_column)))
        temp_name = input_column.tolist()
        new_names = []
        for i in range(len(temp_name)):
            x = temp_dict[temp_name[i]]
            if x > 0:
                new_names.append("_".join([temp_name[i], str(x)]))
            else:
                new_names.append(temp_name[i])
            temp_dict[temp_name[i]] += 1
        single_df.drop("name", axis=1, inplace=True)
        single_df["name"] = new_names


def create_hmm_nodes(all_hmms_df):
    unique_ids = all_hmms_df.copy(deep=True)
    return unique_ids.drop_duplicates("sha512t24u")[["sha512t24u", "model_length"]]


def split_df_by_db(all_hmms_df, db_name_list):
    temp_dict = dict.fromkeys(db_name_list, None)
    # get df indices for each hmm db
    indices_dict = {
        key: all_hmms_df["rel_path"].str.startswith(key)
        for key, value in temp_dict.items()
    }
    # split df by source db
    return {
        key: all_hmms_df.copy(deep=True).loc[value]
        for key, value in indices_dict.items()
    }


def remove_columns_from_dict_of_dfs(dict_of_dfs, cols_to_keep):
    for value in dict_of_dfs.values():
        cols_to_drop = [x for x in list(value.columns.values) if x not in cols_to_keep]
        value.drop(cols_to_drop, axis=1, inplace=True)


def main():
    args = parser.parse_args()
    all_hmms_path = Path(args.all_hmms_path)
    all_hmms_df = pd.read_csv(all_hmms_path, sep="\t")
    # Create a table for the neo4j hmm "nodes" by dropping all rows with duplicate "sha512t24u"
    write_tsv(create_hmm_nodes(all_hmms_df=all_hmms_df), "sg_hmm_nodes_out")

    database_top_dirs = HMM_SOURCES
    # split all_hmm df into dict of dbs
    dict_of_dfs = split_df_by_db(
        all_hmms_df=all_hmms_df, db_name_list=database_top_dirs
    )
    if any([value.empty for value in dict_of_dfs.values()]):
        log.warning("Using subset of available hmms")  # TODO: more helpful message
    del all_hmms_df
    # subset df columns
    remove_columns_from_dict_of_dfs(
        dict_of_dfs=dict_of_dfs,
        cols_to_keep=["sha512t24u", "acc", "name", "description", "category"],
    )

    unique_hmm_names(dict_of_dfs=dict_of_dfs)

    # need to make sure columns are in specific order for import in neo4j
    for key, value in dict_of_dfs.items():
        write_tsv(
            value,
            "".join([key, "_hmms_out"]),
            columns=["sha512t24u", "acc", "name", "description", "category"],
        )


if __name__ == "__main__":
    main()
