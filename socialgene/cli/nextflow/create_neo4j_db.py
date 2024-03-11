import argparse

from socialgene.neo4j.utils.admin_import import Neo4jAdminImport

parser = argparse.ArgumentParser(
    description="Write the header files for neo4j admin import"
)

parser.add_argument(
    "--neo4j_top_dir",
    help="Output directory where header files will be written to",
    required=True,
)

parser.add_argument(
    "--sg_modules",
    metavar="socialgene module list",
    help="List of relationship labels to not create arguments for; e.g. to add mmseqs: base,mmseqs2 ",
    required=True,
    default="base",
    nargs="+",
)

parser.add_argument(
    "--cpus",
    help="Number of cpus to give to neo4j admin import",
    required=False,
    default=4,
)

parser.add_argument(
    "--docker_version",
    help="Sets 'chasemc2/sgnf-sgpy:{docker_version}'",
    required=False,
    default="latest",
)

parser.add_argument(
    "--additional_args",
    help="Additional arguments to pass to the neo4j admin import docker command, given as single string",
    required=False,
    default="",
)
parser.add_argument(
    "--uid",
    type=int,
    help="uid to give to docker, default is the current user",
    required=False,
)
parser.add_argument(
    "--gid",
    type=int,
    help="gid to give to docker, default is the current user",
    required=False,
)
parser.add_argument(
    "--dryrun",
    help="just write to file instead of running",
    required=False,
    default=None,
    action=argparse.BooleanOptionalAction,
)
parser.add_argument(
    "--dryrun_filepath",
    help="Filepath of neo4j import script that will be created",
    required=False,
    default="./build_db.sh",
)
parser.add_argument(
    "--docker",
    help="Use this flag to to add docker-specific commands",
    required=False,
    default=False,
    action=argparse.BooleanOptionalAction,
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
