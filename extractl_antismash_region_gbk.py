#!/usr/bin/env python
# -*- encoding: utf-8 -*-

import argparse
import gzip
import io
import re
import tarfile
from multiprocessing import Pool
from pathlib import Path

from Bio import SeqIO

parser = argparse.ArgumentParser(description="Extract from antismash json")
parser.add_argument(
    "--input_dir",
    metavar="filepath",
    help="path containing many tgz of antismash results",
    required=True,
)

parser.add_argument(
    "--ncpus",
    metavar="int",
    help="ncpus",
    required=True,
)

re_region_gbk = re.compile("region[0-9]+\\.gbk")
re_genomic_gbk = re.compile("genomic.gbk")

OUTDIR = "/home/chase/Downloads/tt/aaa"


def read_write(input_path):
    dbxrefs = []
    assembly_name = "oops"
    # this is complicated because have to add assembly accession back in because it's missing in the region files
    with tarfile.open(input_path, "r") as tar:
        for member in tar:
            if re_genomic_gbk.search(member.name):
                assembly_name = Path(member.name).stem.removesuffix("_genomic")
                with io.TextIOWrapper(tar.extractfile(member)) as handle:
                    for i in SeqIO.parse(handle, "genbank"):
                        dbxrefs = i.dbxrefs
                        continue
        for member in tar:
            if re_region_gbk.search(member.name):
                with gzip.open(
                    Path(OUTDIR, f"{assembly_name}.gbk.gz"), "at"
                ) as output_handle:
                    with io.TextIOWrapper(tar.extractfile(member)) as handle:
                        for record in SeqIO.parse(handle, "genbank"):
                            record.dbxrefs = dbxrefs
                            SeqIO.write(record, output_handle, "genbank")


def find_tgz(input_dir):
    return Path(input_dir).glob("*.tgz")


def kinda(one_path):
    read_write(one_path)


def main():
    args = parser.parse_args()
    # for i in find_tgz(args.input_dir):
    #     read_write(i)

    with Pool(int(args.ncpus)) as pool:
        for i in pool.imap(kinda, find_tgz(args.input_dir)):
            pass


if __name__ == "__main__":
    main()
