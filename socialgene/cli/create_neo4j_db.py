# python dependencies
import argparse

# external dependencies

# internal dependencies
from socialgene.neo4j.admin_import import Neo4jAdminImport
from socialgene.neo4j.sg_modules import parse_hmmlist_input


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
    "--hmmlist",
    metavar="comma sep list",
    help="List of HMM models used",
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

    parsed_hmmlist = parse_hmmlist_input(args.hmmlist)

    temp = Neo4jAdminImport(
        neo4j_top_dir=args.neo4j_top_dir,
        input_sg_modules=args.sg_modules,
        hmmlist=parsed_hmmlist,
        cpus=int(args.cpus),
        additional_args=args.additional_args,
        uid=uid,
        gid=gid,
    )
    if args.dryrun:
        temp.arg_builder(temp.input_sg_modules, temp.hmmlist)
        temp._check_files()
        temp._escape_arg_glob()
        with open("build_db.sh", "w") as handle:
            handle.writelines(" ".join(temp._neo4j_admin_import_args()))
    else:
        temp.run_neo4j_admin_import()


if __name__ == "__main__":
    main()
