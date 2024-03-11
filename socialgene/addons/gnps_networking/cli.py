import argparse
import sys

from socialgene.addons.gnps_networking.parse import GNPS_SNETS
from socialgene.utils.logging import log


def create_arg_parser():
    parser = argparse.ArgumentParser(
        description="Integrate a NPAtlas molecular network into a SocialGene Neo4j Database"
    )
    parser.add_argument(
        "--gnps_dirpath",
        help="Local path to unzipped GNPS Molecular Networking download",
        default=None,
        required=False,
    )

    parser.add_argument(
        "--params_xml_path",
        help="Local path to params.xml file",
        default=None,
        required=False,
    )

    parser.add_argument(
        "--specnets_path",
        help="Local path to specnets file",
        default=None,
        required=False,
    )

    parser.add_argument(
        "--selfloop_path",
        help="Local path to selfloop file",
        default=None,
        required=False,
    )

    parser.add_argument(
        "--clustersummary_path",
        help="Local path to clustersummary file",
        default=None,
        required=False,
    )

    parser.add_argument(
        "--clusterinfo_path",
        help="Local path to clusterinfo file",
        default=None,
        required=False,
    )

    parser.add_argument(
        "--map_path",
        help="Local path to map file (tab separated file with two columns and no header; the first column is the assembly id and the second column is the GNPS mass spec filename)",
        default=None,
        required=False,
    )

    parser.add_argument(
        "--force",
        help="Force integration of GNPS molecular network into the database",
        action="store_true",
        default=False,
        required=False,
    )
    return parser


def main():
    parser = create_arg_parser()
    args = parser.parse_args()
    if args.gnps_dirpath:
        gnps = GNPS_SNETS(gnps_dirpath=args.gnps_dirpath, map_path=args.map_path)
        gnps._find_gnps_paths()
    else:
        # check args must contain all the paths to the GNPS files and they must exist
        for i in [
            args.params_xml_path,
            args.specnets_path,
            args.selfloop_path,
            args.clustersummary_path,
            args.clusterinfo_path,
            # args.map_path,
        ]:
            if not i:
                log.error(
                    "If --gnps_dirpath is not provided, then all GNPS file paths must be provided"
                )
                return
        gnps = GNPS_SNETS(
            params_xml_path=args.params_xml_path,
            specnets_path=args.specnets_path,
            selfloop_path=args.selfloop_path,
            clustersummary_path=args.clustersummary_path,
            clusterinfo_path=args.clusterinfo_path,
            map_path=args.map_path,
        )
    gnps._parse()
    if not args.force:
        if args.map_path:
            # assembly_check = gnps._validate_map()
            # if assembly_check.shape[0] > 0:
            #     if not args.force:
            #         raise ValueError(
            #             "Not all assemblies in the GNPS network were found in the database or vice-versa.\nTo continue anyways, use the --force flag."
            #         )
            ...
        else:
            log.error(
                "To link MS data to genomes, provide --map_path; to skip this step, use the --force flag."
            )
            sys.exit(1)

    # loop through all the methods
    for i in [
        "_add_input_spectra",
        "_create_library_hit_nodes",
        "_add_library_hit_nodes",
        "_add_library_hit_nodes",
        "_add_cluster_nodes",
        "_add_input_ms_files",
        "_link_cluster_to_spectrum",
        "_link_spectrum_to_ms_file",
        "_link_npclassifiers",
        "_link_cluster_to_library_no_spec_id",
        "_link_cluster_to_library",
        "_link_network",
        "_link_library_to_chem",
    ]:
        try:
            getattr(gnps, i)()
        except Exception as e:
            log.error(f"Failed to {i}: {e}")

    if args.map_path:
        try:
            gnps._link_ms_file_to_assembly()
        except Exception as e:
            log.error(f"Failed to link mass spec file to assembly: {e}")

    # gnps.add_input_spectra_to_gnps_clusters()
    # gnps.add_links_between_gnps_clusters()
    log.info(
        "GNPS molecular network has been integrated into the SocialGene Neo4j database"
    )


if __name__ == "__main__":
    main()
