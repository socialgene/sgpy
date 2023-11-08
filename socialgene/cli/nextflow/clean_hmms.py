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
    "--glob",
    default="**/*.socialgene.hmm.gz",
    help="Glob for hmm files",
    type=str,
    required=False,
)


def run_nf_workflow(input_dir, outdir, input_glob="**/*.socialgene.hmm.gz"):
    # glob determined by filepath assigned in the nextflow workflow: https://github.com/socialgene/sgnf/blob/main/bin/hmmconvert_loop.sh
    hmms_object = HmmModelHandler()
    # expects all HMM databases to be saved in a folder named for
    # the database in which it was downloaded from (e.g. antismash, pfam, etc.)
    # so that the name of the source can be assigned
    dirpaths = [
        os.path.join(input_dir, i) for i in os.listdir(input_dir) if i in HMM_SOURCES
    ]
    dirpaths.sort()
    hmm_paths = []
    for dirpath in dirpaths:
        # hmm_paths is a list of dicts: [{filepath:"", base_dir:""}]
        hmm_paths.extend(
            [
                {"filepath": i, "base_dir": dirpath}
                for i in sorted(Path(dirpath).glob(input_glob))
            ]
        )
    # sort for reproducibility
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
        outdir=outdir,
        hash_as_name=True,
    )
    hmms_object.write_metadata_tsv(outdir=outdir, header=False)
    hmms_object.write_hmm_node_tsv(outdir=outdir)


def main():
    args = parser.parse_args()
    input_dir = args.input_dir
    outdir = args.outdir

    run_nf_workflow(
        input_dir=input_dir,
        outdir=outdir,
        input_glob=args.glob,
    )


if __name__ == "__main__":
    main()
