import csv
from pathlib import Path
from typing import List

from socialgene.neo4j.schema.modules import Modules
from socialgene.utils.logging import log


class SocialgeneModules(Modules):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.selected_nodes = set()
        self.selected_relationships = set()

    def add_modules(self, module_list: List[str]):
        """Take a list of modules and create sets of corresponding node/relationship objects

        Args:
            module_list (List[str]): list of socialgene modules (e.g. ["base", "protein"])
        """
        # if only a single module str is provided, make sure it's a list
        for module in list(module_list):
            self.add_module(module)

    def add_module(self, module_id: str):
        if module_id not in self.modules:
            raise ValueError(
                f"Module {module_id} not found. Please select from: {list(self.modules.keys())}"
            )
        else:
            log.debug(f"Adding module: {module_id}")
            for node in self.modules.get(module_id).nodes:
                self.selected_nodes.add(node)
            for rel in self.modules.get(module_id).relationships:
                self.selected_relationships.add(rel)

    def _writer(self, outdir, header, header_filename):
        outpath = Path(outdir, header_filename)
        with open(outpath, "w") as tsv_output_con:
            tsv_writer = csv.writer(
                tsv_output_con,
                delimiter="\t",
                quotechar='"',
                quoting=csv.QUOTE_MINIMAL,
            )
            tsv_writer.writerow(header)
            log.info(f"\tWriting {header_filename} to: {outpath}")

    def write_neo4j_headers(self, outdir: str):
        for node in self.selected_nodes:
            self._writer(
                outdir,
                header=node().header,
                header_filename=node().header_filename,
            )
        for rel in self.selected_relationships:
            self._writer(
                outdir,
                header=rel.header,
                header_filename=rel.header_filename,
            )
