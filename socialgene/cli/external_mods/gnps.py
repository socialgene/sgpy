import argparse
from pathlib import Path

from socialgene.addons.gnps_networking.nr import GNPS_SNETS

parser = argparse.ArgumentParser(
    description="Integrate a GNPS molecular network into a SocialGene Neo4j Database"
)
parser.add_argument(
    "--inputdir",
    help="Path to unzipped directory of GNPS molecular network download",
    default=False,
    required=True,
)


def main():
    args = parser.parse_args()
    if not Path(args.inputdir).exists():
        raise FileNotFoundError(f"Path {args.inputdir} does not exist")
    if not Path(args.inputdir).is_dir():
        raise ValueError(f"Path {args.inputdir} is not a directory")
    gnps = GNPS_SNETS(gnps_dirpath=args.inputdir)
    gnps.add_gnps_library_spectra()
    gnps.add_gnps_library_spectrum_classifications()
    gnps.add_gnps_library_spectrum_ion_sources()
    gnps.add_gnps_library_spectrum_instruments()
    gnps.add_gnps_library_spectrum_organisms()
    gnps.add_gnps_clusters()
    gnps.add_links_between_gnps_clusters_and_assemblies()
    gnps.add_links_between_gnps_clusters()


if __name__ == "__main__":
    main()
