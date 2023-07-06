#
# import argparse

#

#
# from socialgene.utils.logging import log
# from socialgene.parsers.ncbi_feature_table import NCBIFeatureTable
# from socialgene.config import env_vars

# parser = argparse.ArgumentParser(
#     description="Export header files for Neo4j admin import"
# )

# parser.add_argument(
#     "--input",
#     metavar="filepath",
#     help="Input feature tables",
#     required=True,
# )

# parser.add_argument(
#     "--sqlite_path",
#     metavar="filepath",
#     help="Path to SQLite database",
#     required=True,
# )
# parser.add_argument(
#     "--outdir",
#     metavar="filepath",
#     help="Output directory filepath",
#     required=True,
# )

# input_path = (
#     "/home/chase/Documents/socialgene/test_data/BGC0000001_feature_table.txt.gz"
# )


# def main():
#     args = parser.parse_args()
#     outdir = args.outdir
#     log.info(f"Socialgene variables: \n{env_vars}")
#     ft_object = NCBIFeatureTable()
#     ft_object.connect("/home/chase/Documents/socialgene/test_data/hashid.sqlite")
#     ft_object.parse(args.input)

#     ft_object.export_locus_to_protein(outdir=outdir)
#     ft_object.export_protein_info(outdir=outdir)
#     ft_object.export_assembly_to_locus(outdir=outdir)
#     ft_object.export_loci(outdir=outdir)
#     ft_object.export_assembly(outdir=outdir)


# if __name__ == "__main__":
#     main()
