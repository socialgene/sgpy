import pandas as pd


def write_tsv(input_data, filepath, columns=None, mode="w"):
    input_data.to_csv(
        path_or_buf=filepath,
        sep="\t",
        na_rep="",
        float_format=None,
        columns=columns,
        header=False,
        index=False,
        index_label=None,
        mode=mode,
        encoding=None,
        compression="infer",
        quoting=None,
        quotechar='"',
        lineterminator=None,
        chunksize=None,
        date_format=None,
        doublequote=True,
        escapechar=None,
        decimal=".",
        errors="strict",
        storage_options=None,
    )


def read_tsv(filepath, header=None):
    return pd.read_csv(filepath, sep="\t", header=header)
