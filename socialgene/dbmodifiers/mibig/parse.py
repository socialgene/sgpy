import json
import tarfile
import tempfile
from pathlib import Path

import requests
from rich.progress import Progress, SpinnerColumn

from socialgene.dbmodifiers.mibig.compound import Mibig_Compound
from socialgene.neo4j.neo4j import GraphDriver

MIBIG_URL = "https://dl.secondarymetabolites.org/mibig/mibig_json_3.1.tar.gz"


def dl_json():
    with tempfile.TemporaryDirectory() as tmpdirname:
        with open(Path(tmpdirname, "mibig_json_3.1.tar.gz"), "wb") as out_file:
            content = requests.get(MIBIG_URL, stream=True).content
            _ = out_file.write(content)

        tar_path = Path(tmpdirname, "mibig_json_3.1.tar.gz")

        tarfile_object = tarfile.open(tar_path)
        bgc_files = tarfile_object.getnames()
        for i in bgc_files:
           if i.endswith(".json"):
               yield i


single_file = "/home/chase/Downloads/mibig_json_3.1/BGC0001850.json"


with open(single_file, "r") as f:
    data = json.load(f)





class MibigParse:
    def __init__(self, entry):
        self.entry = entry

    def get_class(self):
        try:
            for i in self.entry['cluster']['biosyn_class']:
                self.mibig_biosynthetic_class = Mibig_Biosynthetic_Class(properties={"uid": i})
        except Exception as e:
            print(e)

    def get_compounds(self):
        try:
            for i in self.entry['cluster']['compounds']:
                self.mibig_compound = Mibig_Compound(properties={"uid": i})
        except Exception as e:
            print(e)


for i in data['cluster']['biosyn_class']:
    Mibig_Biosynthetic_Class(properties={"uid": i})



data['cluster']['compounds'][0]



















if __name__ == "__main__":
    add_mibig_info_to_neo4j()

