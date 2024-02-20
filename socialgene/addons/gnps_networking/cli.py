
import argparse

from socialgene.addons.gnps_networking.parse import GNPS_SNETS
from socialgene.utils.logging import log


def create_arg_parser():
    """ "Creates and returns the ArgumentParser object."""
    parser = argparse.ArgumentParser(
        description="Integrate a NPAtlas molecular network into a SocialGene Neo4j Database"
    )
    parser.add_argument(
        "--input",
        help="Local path to unzipped GNPS Molecular Networking download",
        default=None,
        required=True,
    )
    parser.add_argument(
        "--regex",
        help="Pattern of MS file names. If '_' then it expects 'assembly-id_*.mzML'; if '__' then 'assembly-id__*.mzML'; can also be a custom regex pattern. Defaults to '_'.",
        default="_",
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
    input = args.input
    gnps = GNPS_SNETS(gnps_dirpath=input)
    gnps._filename_to_assembly(regex=args.regex)
    assembly_check = gnps._check_db_for_assemblies()
    if not assembly_check:
        if not args.force:
            log.warning("Not all assemblies in the GNPS network were found in the database or vice-versa.\nTo continue anyways, use the --force flag.")
            return
    gnps._library_hit_nodes()
    gnps.add_library_hit_nodes_to_neo4j()
    gnps.link_library_to_chem()
    gnps.create_cluster_nodes(create=False)
    gnps.link_npclassifiers()
    gnps.link_cluster_to_library()
    gnps.link_network()
    gnps.link_assembly_to_cluster(create=False)
    # gnps.add_input_spectra_to_gnps_clusters()
    # gnps.add_links_between_gnps_clusters()
    log.info("GNPS molecular network has been integrated into the SocialGene Neo4j database")
