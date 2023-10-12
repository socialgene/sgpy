import json
import tarfile
import tempfile
from pathlib import Path

import requests
from rich.progress import Progress, SpinnerColumn

from socialgene.dbmodifiers.mibig.compound import Mibig_Compound
from socialgene.neo4j.neo4j import GraphDriver

# #


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
        with Progress(
            SpinnerColumn(spinner_name="runner"),
            *Progress.get_default_columns(),
        ) as progress:
            task = progress.add_task(
                "Updating MIBiG nodes and relationships...", total=len(bgc_files)
            )
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
                progress.update(task, advance=1)


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
