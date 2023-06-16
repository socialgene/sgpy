from typing import List

import csv
from pathlib import Path

from socialgene.neo4j.schema.define_modules import Modules
from socialgene.neo4j.schema.define_nodes import Nodes
from socialgene.neo4j.schema.define_relationships import Relationships
from socialgene.utils.logging import log


class SocialgeneModules(Modules):
    def __init__(
        self,
    ):
        super().__init__()
        self.nodes = set()
        self.relationships = set()
        self.all_nodes = Nodes()
        self.all_relationships = Relationships()

    def add_modules(self, module_list: List[str]):
        """Take a list of modules and add the corresponding nodes/relationships

        Args:
            module_list (List[str]): list of socialgene modules (e.g. ["base", "protein"])
        """
        if isinstance(module_list, str):
            module_list = [module_list]

        for module in module_list:
            if module in self.modules:
                module_def = self.modules.get(module)
                for node in self.all_nodes.get_nodes(module_def.nodes):
                    self.nodes.add(node)
                for rel in list(
                    self.all_relationships.get_relationships(module_def.relationships)
                ):
                    self.relationships.add(rel)

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
        for node in self.nodes:
            self._writer(
                outdir, header=node.header, header_filename=node.header_filename
            )
        for rel in self.relationships:
            self._writer(outdir, header=rel.header, header_filename=rel.header_filename)
