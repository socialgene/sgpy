import argparse
import os
from pathlib import Path

from rich.progress import Progress, SpinnerColumn

from socialgene.parsers.hmmmodel import HMM_SOURCES, HmmModelHandler
from socialgene.utils.logging import log

parser = argparse.ArgumentParser(description="Parse NcbiAssembliessdsd taxonomy")
parser.add_argument(
    "--input_dir",
    metavar="input filepath",
    type=Path,
    help="input_dir",
    required=True,
)
parser.add_argument(
    "--outdir",
    metavar="output directory",
    type=Path,
    help="outdir",
    required=True,
)
parser.add_argument(
    "--numoutfiles",
    metavar="Number of out files",
    type=int,
    help="numoutfiles file",
    required=True,
)
parser.add_argument(
    "--splitcutoffs",
    default="false",
    help="Create two file, one with models that have cutoff values, one with models that don't",
    required=True,
    action=argparse.BooleanOptionalAction,
)


def run_nf_workflow(
    input_dir, outdir, n_files, splitcutoffs, input_glob="**/*.socialgene.hmm.gz"
):
    # glob determined by https://github.com/socialgene/sgnf/blob/17f9cde454c13a91bac95917237cf35ddcbe52c3/bin/hmmconvert_loop.shL10
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
        transient=True,
    ) as progress:
        task = progress.add_task("Progress...", total=len(hmm_paths))
        for i in hmm_paths:
            hmms_object.read(filepath=i.get("filepath"), base_dir=i.get("base_dir"))
            log.debug(f"Parsing {i.get('filepath')}\r")
            progress.update(task, advance=1)
    # Update info, make non-redundant, and write out
    hmms_object.hydrate_cull()
    hmms_object.remove_duplicate_and_old_pfam()
    hmms_object.remove_duplicate_hash()
    hmms_object.write(
        outdir=outdir, n_files=n_files, hash_as_name=True, splitcutoffs=splitcutoffs
    )
    hmms_object.write_metadata_tsv(outdir=outdir, header=False)
    hmms_object.write_hmm_node_tsv(outdir=outdir)


def main():
    args = parser.parse_args()
    input_dir = args.input_dir
    outdir = args.outdir
    numoutfiles = int(args.numoutfiles)

    run_nf_workflow(
        input_dir=input_dir,
        outdir=outdir,
        n_files=numoutfiles,
        splitcutoffs=args.splitcutoffs,
    )


if __name__ == "__main__":
    main()
