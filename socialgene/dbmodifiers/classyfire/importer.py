import argparse

from rich.progress import Progress, SpinnerColumn

from socialgene.external_db_classes.classyfire import ClassyFire
from socialgene.utils.logging import log

parser = argparse.ArgumentParser(
    description="Download and add ClassyFire ontology to a running SocialGene database"
)


def main():
    obs = ClassyFire.download()
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
