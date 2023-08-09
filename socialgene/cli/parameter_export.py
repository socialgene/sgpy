import argparse
import csv
import platform

from socialgene.config import env_vars
from socialgene.utils.logging import log

parser = argparse.ArgumentParser(
    description="Export common parameters, especially for nextflow/neo4j"
)

parser.add_argument(
    "--outpath",
    help="Output filepath",
    default="socialgene_parameters",
)

# expect a string "None"; not None
parser.add_argument("--genome_download_command", required=False, default="None")

fields = {
    "SG_LOC_NEO4J": None,
    "SG_LOC_HMMS": None,
    "NEO4J_dbms_memory_pagecache_size": None,
    "NEO4J_dbms_memory_heap_initial__size": None,
    "NEO4J_dbms_memory_heap_max__size": None,
    "HMMSEARCH_IEVALUE": None,
    "HMMSEARCH_BACKGROUND": None,
    "HMMSEARCH_BIASFILTER": None,
    "HMMSEARCH_NULL2": None,
    "HMMSEARCH_SEED": None,
    "HMMSEARCH_Z": None,
    "HMMSEARCH_DOMZ": None,
    "HMMSEARCH_F1": None,
    "HMMSEARCH_F2": None,
    "HMMSEARCH_F3": None,
    "HMMSEARCH_E": None,
    "HMMSEARCH_DOME": None,
    "HMMSEARCH_INCE": None,
    "HMMSEARCH_INCDOME": None,
    "HMMSEARCH_BITCUTOFFS": None,
}


def main():
    args = parser.parse_args()

    save_list = ["socialgene_parameter_export"]

    with open(args.outpath, "a") as f:
        for k, v in env_vars.items():
            if k in fields:
                fields[k] = v
            else:
                log.info(f"Unexpected parameter: '{k}'")
        for v in fields.values():
            save_list.append(v)

        # TODO: trycatch all this and check length at end
        save_list.append(platform.sys.platform)
        save_list.append(" ".join(platform.architecture()))
        save_list.append(platform.sys.executable)
        save_list.append(platform.sys.version)
        save_list.append(args.genome_download_command)
        save_list = [str(i) for i in save_list]
        save_list = [i.strip('"') for i in save_list]
        save_list = [i.replace("\n", "") for i in save_list]
        save_list = [i.replace("\t", "") for i in save_list]

        writer = csv.writer(
            f,
            quoting=csv.QUOTE_MINIMAL,
            delimiter="\t",
        )
        writer.writerows([save_list])


# TODO:
# date/time
# Py/Next/Neo versions
# SG version and/or git sha
# domputer info
# source genome info


if __name__ == "__main__":
    main()
