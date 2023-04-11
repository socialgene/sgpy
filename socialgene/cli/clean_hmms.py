# python dependencies
from pathlib import Path
import os
import argparse

# external dependencies
from rich.progress import Progress, SpinnerColumn

# internal dependencies
from socialgene.parsers.hmmmodel import HMM_SOURCES, HmmModelHandler


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
    input_glob = "**/*_socialgenehmm.gz"
    hmms_object = HmmModelHandler()
    # currently expects all HMM databases to be saved in a folder named for
    # the database in which it was downloaded from. To add more, see HMMParser().hmm_dbs
    # Because information about HMM function is sometimes held in directory
    # structure (e.g. antismash)
    dirpaths = [
        os.path.join(input_dir, i) for i in os.listdir(input_dir) if i in HMM_SOURCES
    ]
    hmm_paths = []
    for dirpath in dirpaths:
        # hmm_paths is a list of dicts: [{filepath:"", base_dir:""}]
        hmm_paths.extend(
            [
                {"filepath": i, "base_dir": dirpath}
                for i in Path(dirpath).glob(input_glob)
            ]
        )
    # read and parse all models
    with Progress(
        SpinnerColumn(spinner_name="runner"),
        *Progress.get_default_columns(),
    ) as progress:
        task = progress.add_task("Progress...", total=len(hmm_paths))
        for i in hmm_paths:
            print(i)
            hmms_object.read(filepath=i.get("filepath"), base_dir=i.get("base_dir"))
            progress.console.print(f"Parsing {i.get('filepath')}\r")
            progress.update(task, advance=1)
    # Update info, make non-redundant, and write out
    hmms_object.hydrate_cull()
    hmms_object.remove_duplicate_and_old_pfam()
    hmms_object.remove_duplicate_hash()
    hmms_object.write_culled(outdir=outdir, n_files=n_files, hash_as_name=True)
    hmms_object.write_metadata_tsv(outdir=outdir, header=False)
    hmms_object.write_hmm_node_tsv(outdir=outdir)


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
