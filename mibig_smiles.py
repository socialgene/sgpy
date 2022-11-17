# import json
# from pathlib import Path
# from rich import inspect


# json_dir = "/home/chase/Downloads/mibig_json_3.1"

# pathlist = Path(json_dir).glob("BGC*.json")


# collect_dict = {}

# for json_path in pathlist:
#     try:
#         with open(json_path, "r") as f:
#             data = json.load(f)
#         if "cluster" in data:
#             if "compounds" in data["cluster"]:
#                 # only get the first one
#                 collect_dict[data["cluster"]["mibig_accession"]] = data["cluster"][
#                     "compounds"
#                 ][0]["chem_struct"]
#     except:
#         pass


# from socialgene.neo4j.neo4j import GraphDriver  # grab the the neo4j connection
# import logging

# logging.getLogger("neo4j").setLevel(logging.WARNING)
# logging.getLogger().setLevel(logging.INFO)

# query = """
# WITH $params as inputs
# UNWIND inputs AS input
# MATCH (a1:assembly)
# WHERE a1.id = input[0]
# SET a1.SMILES = input[1]
# """

# with GraphDriver() as db:
#     results = db.run(query, params=[[k, v] for k, v in collect_dict.items()])
