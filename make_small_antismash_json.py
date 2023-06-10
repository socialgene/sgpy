#!/usr/bin/env python
# -*- encoding: utf-8 -*-
import argparse
import json
import tarfile
from multiprocessing import Pool
from pathlib import Path

#  find $myf -name '*.tgz' -exec tar  --to-stdout --wildcards "*_genomic.json" -xf {} \; | jq --compact-output '[.records[].dbxrefs][0] as $db | [{(.records[].id): .records[].areas}] as $loci | {dbxrefs: $db, records: $loci}' > all_antismash.json
# jq --compact-output '[.records[].dbxrefs][0] as $db | [{(.records[].id): .records[].areas}] as $loci | {dbxrefs: $db, records: $loci}'  $myf


parser = argparse.ArgumentParser(description="Extract from antismash json")
parser.add_argument(
    "--input_dir",
    metavar="filepath",
    help="path containing many tgz of antismash results",
    required=True,
)
parser.add_argument(
    "--outpath",
    metavar="filepath",
    help="path to write jsonl",
    required=True,
)

parser.add_argument(
    "--ncpus",
    metavar="int",
    help="ncpus",
    required=True,
)


def read_json(input_path):
    with tarfile.open(input_path, "r") as tar:
        for member in tar:
            if member.name.endswith("_genomic.json"):
                f = tar.extractfile(member)
                return json.load(f).get("records")


def find_tgz(input_dir):
    return Path(input_dir).glob("*.tgz")


def bro(records):
    assembly = None
    try:
        for i in records[0].get("dbxrefs"):
            if i.startswith("Assembly:"):
                assembly = i.replace("Assembly:", "")
    except Exception:
        pass
    return {
        "assembly": assembly,
        "records": [
            {
                record.get("id"): record.get("areas")
                for record in records
                if record.get("areas")
            }
        ],
    }


def kinda(one_path):
    return bro(read_json(one_path))


def main():
    args = parser.parse_args()
    with Pool(int(args.ncpus)) as pool:
        with open(args.outpath, "w") as outfile:
            for result in pool.imap(kinda, find_tgz(args.input_dir)):
                json.dump(result, outfile)
                outfile.write("\n")


if __name__ == "__main__":
    main()
