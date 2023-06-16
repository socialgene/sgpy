#
# from urllib.request import urlopen
# from pathlib import Path

#
# import xml.etree.ElementTree as ET

#


# # these are run manually to write lines in a file that are the attributes present
# # in NCBI's biosample API xml return e.g.

# # tree = read_biosample_attributes_xml()
# # attribute_list = parse_biosample_attributes_xml(tree=tree)
# # write_attributes(
# #     outdir="/home/chase/Documents/socialgene/socialgene/src/socialgene/data",
# #     attribute_list=attribute_list,
# # )


# def read_biosample_attributes_xml():

#     # read the xml format from ncbi
#     var_url = urlopen(
#         "https://www.ncbi.nlm.nih.gov/biosample/docs/attributes/?format=xml"
#     )
#     return ET.parse(var_url)


# def parse_biosample_attributes_xml(tree):
#     """Parse the returned XML

#     Args:
#         tree (xml.etree.ElementTree.ElementTree): result of read_biosample_attributes_xml()
#     Returns:
#         list:
#     """

#     root = tree.getroot()
#     attr = []
#     for i in root.iterfind("Attribute/HarmonizedName"):
#         attr.append(i.text)
#     # standardize output
#     attr = set(attr)
#     attr = list(attr)
#     attr.sort()
#     return attr


# def write_attributes(outdir, attribute_list):
#     """Write the attributes to a file

#     Args:
#         outdir (str): directory where a file named "biosample_attributes" will be written
#         attribute_list (list): attributes to write to the file
#     """
#     biosample_attribute_tags_path = Path(outdir, "biosample_attributes")
#     with open(biosample_attribute_tags_path, "w") as handle:
#         handle.write("\n".join(attribute_list))
