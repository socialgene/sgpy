import json
from pathlib import Path
import argparse
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

parser = argparse.ArgumentParser(description="Add npatlas info to neo4j database")
parser.add_argument(
    "--input_dir",
    metavar="filepath",
    help="input_dir",
    required=True,
)

URL = "https://dl.secondarymetabolites.org/mibig/mibig_json_3.1.tar.gz"



npatlas_path = "/home/chase/Downloads/NPAtlas_download.json"
import json

with open(npatlas_path, "r") as h:
    a =json.load(h)


from socialgene.dbmodifiers.npatlas.npatlas import Npatlas


z=[Npatlas(i) for i in a]

from rich import inspect
inspect([i for i in z if i.uid == 'NPA025463'][0])
[i for i in z if i.uid == 'NPA025463'][0].gnps


inspect([i for i in z if i.mibig][0])



    external_descriptors
        "source": "CHEBI",
          "source": "KEGG",
          "source": "LIPID MAPS",
          "source": "META CYC",
