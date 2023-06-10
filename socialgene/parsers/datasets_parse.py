#
# from collections import OrderedDict
# from pathlib import Path
# import json
# import csv

#

#
# from socialgene.utils.logging import log


# class DatasetsParse:
#     pass


# class Biosample:

#     biosample_attributes = [
#         "biosample_accession",
#         "collection_date",
#         "culture_collection",
#         "env_broad_scale",
#         "env_local_scale",
#         "env_medium",
#         "env_package",
#         "geo_loc_name",
#         "host",
#         "investigation_type",
#         "isol_growth_condt",
#         "isolation_source",
#         "lat_lon",
#         "num_replicons",
#         "project_name",
#         "ref_biomaterial",
#         "sample_type",
#         "type-material",
#     ]

#     def __init__(self, jsonl_path, outdir):
#         self.jsonl_path = jsonl_path
#         self.outdir = outdir
#         self.jsonl_dict = None

#         self._validate()
#         # self._load_jsonl()
#         self._parse_biosamples()

#     def _validate(self):
#         for i in ["assembly_to_biosample", "biosamples"]:
#             temp_path = Path(self.outdir, i)
#             if temp_path.exists() and temp_path.stat().st_size > 0:
#                 raise IOError(
#                     f"File exists already and isn't empty, taking the coward's way out and aborting: {temp_path}"
#                 )

#     def _get_all_biosample_attribute_names(self):
#         # this may be unnecessary but I don't know if the fields change between assemblies
#         self.biosample_attributes = set()
#         for i in self.jsonl_dict:
#             for att in i["assemblyInfo"]["biosample"]["attributes"]:
#                 self.biosample_attributes.add(att["name"])
#         self.biosample_attributes = list(self.biosample_attributes)
#         self.biosample_attributes.sort()
#         # the results of self.biosample_attributes is hardcoded as a dictionary for ensured consistency

#     def _parse_biosamples(self):
#         # File reading/writing here top prevent having data copies everywhere

#         with open(self.jsonl_path) as input_data:
#             for line in input_data:
#                 jsonl_dict = json.loads(line)
#                 with open(
#                     Path(self.outdir, "assembly_to_biosample"), "a"
#                 ) as assembly_to_biosample_handle, open(
#                     Path(self.outdir, "biosamples"), "a"
#                 ) as biosamples_handle:
#                     assembly_to_biosample_tsv_writer = csv.writer(
#                         assembly_to_biosample_handle,
#                         delimiter="\t",
#                         quotechar='"',
#                         quoting=csv.QUOTE_MINIMAL,
#                     )
#                     biosamples_tsv_writer = csv.writer(
#                         biosamples_handle,
#                         delimiter="\t",
#                         quotechar='"',
#                         quoting=csv.QUOTE_MINIMAL,
#                     )

#                     assembly_to_biosample_tsv_writer.writerow(
#                         [
#                             jsonl_dict["assemblyInfo"]["assemblyAccession"],
#                             jsonl_dict["assemblyInfo"]["biosample"]["accession"],
#                         ]
#                     )
#                     biosample_attribute_dict = OrderedDict(
#                         {i: None for i in self.biosample_attributes}
#                     )
#                     biosample_attribute_dict["biosample_accession"] = jsonl_dict[
#                         "assemblyInfo"
#                     ]["biosample"]["accession"]
#                     for att in jsonl_dict["assemblyInfo"]["biosample"]["attributes"]:
#                         try:
#                             temp = att["value"]
#                             temp.replace("\t", " ")
#                             if not att["value"] == "missing":
#                                 biosample_attribute_dict[att["name"]] = att["value"]
#                             if att["name"] == "culture_collection":
#                                 try:
#                                     temp = att["value"].split
#                                 except Exception as e:
#                                     log.error(e)
#                         except Exception as e:
#                             log.debug(f"{att['name']} not found in dict")
#                             log.error(e)
#                             pass
#                     log.info(list(biosample_attribute_dict.values()))
#                     biosamples_tsv_writer.writerow(
#                         list(biosample_attribute_dict.values())
#                     )
