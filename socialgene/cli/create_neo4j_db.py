# python dependencies
import argparse

# external dependencies

# internal dependencies
from socialgene.neo4j.admin_import import Neo4jAdminImport

from socialgene.parsers.hmmmodel import HMM_SOURCES, parse_hmmlist_input

parser = argparse.ArgumentParser(
    description="Write the header files for neo4j admin import"
)

parser.add_argument(
    "--neo4j_top_dir",
    metavar="filepath",
    help="Output directory where header files will be written to",
    required=True,
)

parser.add_argument(
    "--sg_modules",
    metavar="comma sep list",
    help="List of relationship labels to not create arguments for; e.g. to add mmseqs: base,mmseqs2 ",
    required=True,
    default="base",
    nargs="+",
)

parser.add_argument(
    "--cpus",
    metavar="int",
    help="Number of cpus to give to neo4j admin import",
    required=False,
    default=4,
)

parser.add_argument(
    "--additional_args",
    metavar="str",
    help="Additional arguments to pass to the neo4j admin import docker command, given as single string",
    required=False,
    default="",
)
parser.add_argument(
    "--uid",
    metavar="int",
    help="uid to give to docker, default is the current user",
    required=False,
    default=None,
)
parser.add_argument(
    "--gid",
    metavar="int",
    help="gid to give to docker, default is the current user",
    required=False,
    default=None,
)
parser.add_argument(
    "--dryrun",
    metavar="bool",
    help="just write to file instead of running ",
    required=False,
    default=None,
)
parser.add_argument(
    "--dryrun_filepath",
    metavar="filepath",
    help="Filepath of neo4j import script that will be created",
    required=False,
    default="./build_db.sh",
)
parser.add_argument(
    "--docker",
    metavar="bool",
    help="Use this flag to to add docker-specific commands",
    required=False,
    default=False,
)


def main():
    args = parser.parse_args()

    if args.uid == "None":
        uid = None
    else:
        uid = args.uid

    if args.gid == "None":
        gid = None
    else:
        gid = args.gid

    nai_obj = Neo4jAdminImport(
        neo4j_top_dir=args.neo4j_top_dir,
        module_list=args.sg_modules,
        cpus=int(args.cpus),
        additional_args=args.additional_args,
        uid=uid,
        gid=gid,
    )
    tmp = nai_obj.run(docker=args.docker, dryrun=args.dryrun)
    if args.dryrun:
        with open(args.dryrun_filepath, "w") as handle:
            handle.writelines(" \\\n\t".join(tmp))


if __name__ == "__main__":
    main()
