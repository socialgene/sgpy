#!/usr/bin/env python

# python dependencies
from pathlib import Path
import os
import argparse

# external dependencies
from rich.progress import Progress, SpinnerColumn

# internal dependencies
from socialgene.parsers.hmm import HMMParser


parser = argparse.ArgumentParser(description="Parse NcbiAssembliessdsd taxonomy")
parser.add_argument(
    "--input_dir",
    metavar="filepath",
    help="input_dir",
    required=True,
)
parser.add_argument(
    "--outdir",
    metavar="filepath",
    help="outdir",
    required=True,
)
parser.add_argument(
    "--numoutfiles",
    metavar="int",
    help="numoutfiles file",
    required=True,
)


def run_nf_workflow(input_dir, outdir, n_files):
    hmms_object = HMMParser()
    # currently expects all HMM databases to be saved in a folder named for
    # the database in which it was downloaded from. To add more, see HMMParser().hmm_dbs
    # Because information about HMM function is sometimes held in directory
    # structure (e.g. antismash)
    expected_absolute_dirs = {
        k: os.path.join(input_dir, k) for k in HMMParser().hmm_dbs
    }
    # Change strings to Paths
    expected_absolute_dirs_pathlib = {
        key: Path(value) for (key, value) in expected_absolute_dirs.items()
    }
    # Find the up-converted hmm files
    temp_paths = {
        # instead of using the generator from rglob, sort the list of paths it
        # creates so files will always be read in same order
        key: sorted(value.rglob("**/*_socialgene*"))
        for (key, value) in expected_absolute_dirs_pathlib.items()
    }
    # read and parse all models
    with Progress(
        SpinnerColumn(spinner_name="runner"),
        *Progress.get_default_columns(),
    ) as progress:
        task = progress.add_task(
            "Progress...", total=sum([len(i) for i in temp_paths.values()])
        )
        for filename, models in temp_paths.items():
            for path in models:
                hmms_object.parse_single_model_file(
                    input_dir=input_dir, input_path=path, model_source=filename
                )
                progress.console.print(f"Parsing {filename}: {path.stem}\r")
                progress.update(task, advance=1)
    # Update info, make non-redundant, and write out
    hmms_object.change_model_name_to_hash()
    hmms_object.create_nr_hmm_dict()
    hmms_object.write_hmms(outdir=outdir, n_files=n_files)
    hmms_object.write_all_hmm_info(outdir=outdir)


def main():
    args = parser.parse_args()
    # input_dir = "/home/chase/Documents/socialgene_outdir/hmm_downloads"
    # numoutfiles = 5
    # outdir = "/home/chase/Documents/socialgene_outdir/socialgene_hmms_info"
    input_dir = args.input_dir
    outdir = args.outdir
    numoutfiles = int(args.numoutfiles)

    run_nf_workflow(input_dir=input_dir, outdir=outdir, n_files=numoutfiles)


if __name__ == "__main__":
    main()
