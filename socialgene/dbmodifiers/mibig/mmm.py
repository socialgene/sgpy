import json
from pathlib import Path
import tempfile
import requests
import tarfile
from socialgene.dbmodifiers.mibig.nrps import NRPS
from socialgene.neo4j.neo4j import GraphDriver
from socialgene.utils.logging import log
from socialgene.config import env_vars
from socialgene.dbmodifiers.mibig.compound import Mibig_Compound
from rich import inspect
from neo4j import GraphDatabase
from socialgene.config import env_vars
import logging

# logging.getLogger("neo4j").setLevel(logging.WARNING)
# logging.getLogger().setLevel(logging.INFO)


# CREATE CONSTRAINT FOR (n:chemical) REQUIRE n.uid IS UNIQUE;


MIBIG_URL = "https://dl.secondarymetabolites.org/mibig/mibig_json_3.1.tar.gz"


def add_mibig_info_to_neo4j():
    with tempfile.TemporaryDirectory() as tmpdirname:
        with open(Path(tmpdirname, "mibig_json_3.1.tar.gz"), "wb") as out_file:
            content = requests.get(MIBIG_URL, stream=True).content
            _ = out_file.write(content)

        tar_path = Path(tmpdirname, "mibig_json_3.1.tar.gz")

        tarfile_object = tarfile.open(tar_path)
        bgc_files = tarfile_object.getnames()
        bgc_files = [i for i in bgc_files if i.endswith(".json")]

        for file in bgc_files:
            input_json = json.load(tarfile_object.extractfile(file))
            for cmpd in input_json["cluster"]["compounds"]:
                bro = Mibig_Compound(cmpd)
                bro.assign()
                bro.write_all()
                (
                    GraphDriver().driver.execute_query(
                        """
                        MATCH (a1:assembly {uid: $bgc_id})
                        MATCH (chem:chemical {uid: $chem_uid})
                        MERGE (a1)-[:PRODUCES]->(chem)
                        """,
                        chem_uid=bro.uid,
                        bgc_id=input_json["cluster"]["mibig_accession"],
                        database_="neo4j",
                    )
                )


if __name__ == "__main__":
    add_mibig_info_to_neo4j()

# for file in bgc_files:
#     input_json = json.load(tarfile_object.extractfile(file))
#     if "nrp" in input_json["cluster"]:
#         if "cyclic" in input_json["cluster"]["nrp"]:
#             raise
#         input_json["cluster"]["polyketide"]


# input_json["cluster"]["polyketide"]["synthases"][0]["subclass"]


# for file in bgc_files:
#     input_json = json.load(tarfile_object.extractfile(file))
#     a = NRPS()
#     a.assign(input_json)
#     a.write_to_neo4j()
