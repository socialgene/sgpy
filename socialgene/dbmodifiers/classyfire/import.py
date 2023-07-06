import argparse
import requests
import argparse

import requests
from zipfile import ZipFile
from io import BytesIO
import logging
from socialgene.utils.logging import log

from socialgene.external_db_classes.classyfire import ClassyFire
from rich.progress import Progress, SpinnerColumn

logging.getLogger("neo4j").setLevel(logging.WARNING)

parser = argparse.ArgumentParser(
    description="Download and add ClassyFire ontology to a running SocialGene database"
)

OBO_URL = (
    "http://classyfire.wishartlab.com/system/downloads/1_0/chemont/ChemOnt_2_1.obo.zip"
)


def main():
    response = requests.get(OBO_URL)
    response = requests.get(OBO_URL, stream=True)
    obs = []
    log.info(f"Downloading and parsing data from:\n\t{OBO_URL}")
    with ZipFile(BytesIO(response.content)).open("ChemOnt_2_1.obo") as h:
        term_open = False
        for line in h:
            line = line.decode("utf-8")
            if line.startswith("[Term]"):
                term_open = True
                obj = ClassyFire()
            if not term_open:
                continue
            obj.assign(line)
            if line == "\n":
                term_open = False
                obs.append(obj)

    log.info(f"Modifying Neo4j database with {len(obs)} terms")

    with Progress(
        SpinnerColumn(spinner_name="runner"),
        *Progress.get_default_columns(),
    ) as progress:
        task = progress.add_task(
            "Creating CHEMONT and related CHEBI nodes...", total=len(obs)
        )
        for obj in obs:
            obj._create_chemont_nodes()
            obj.create_chebi_nodes()
            progress.update(task, advance=1)

    with Progress(
        SpinnerColumn(spinner_name="runner"),
        *Progress.get_default_columns(),
    ) as progress:
        task = progress.add_task(
            "Linking CHEMONT-to-CHEMONT and CHEMONT-to-CHEBI nodes...", total=len(obs)
        )
        for obj in obs:
            obj.connect_chemont_chebi_nodes()
            obj.connect_chemont_is_a_nodes()
            progress.update(task, advance=1)


if __name__ == "__main__":
    main()
