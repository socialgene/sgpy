#
# from urllib.request import urlopen
# from pathlib import Path

#
# import xml.etree.ElementTree as ET

#


# def read_biosample_attributes_xml():

#     # read the xml format from ncbi
#     var_url = urlopen(
#         "https://www.ncbi.nlm.nih.gov/biosample/docs/attributes/?format=xml"
#     )
#     return ET.parse(var_url)


# def parse_biosample_attributes_xml(tree):
#     """Parse the returned XML

#     Args:
#         tree (xml.etree.ElementTree.ElementTree): [description]

#     Returns:
#         list:
#     """

#     root = tree.getroot()
#     attr = []
#     for i in root.iterfind("Attribute/HarmonizedName"):
#         attr.append(i.text)
#     attr = set(attr)  # make unique
#     attr = list(attr)
#     attr.sort()
#     return attr


# def write_attributes(outdir, attribute_list):
#     biosample_attribute_tags_path = Path(outdir, "biosample_attributes")
#     with open(biosample_attribute_tags_path, "w") as handle:
#         handle.write("\n".join(attribute_list))


# tree = read_biosample_attributes_xml()

# attribute_list = parse_biosample_attributes_xml(tree=tree)

# write_attributes(
#     outdir="/home/chase/Documents/socialgene/socialgene/src/socialgene/data",
#     attribute_list=attribute_list,
# )
