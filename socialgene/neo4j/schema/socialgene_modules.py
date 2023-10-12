import csv
from pathlib import Path
from typing import List

from socialgene.neo4j.schema.modules import Modules
from socialgene.neo4j.schema.nodes import Nodes
from socialgene.neo4j.schema.relationships import Relationships
from socialgene.utils.logging import log


class SocialgeneModules(Modules, Nodes, Relationships):
    def __init__(self, **kwargs):
        super().__init__()
        self.selected_nodes = set()
        self.selected_relationships = set()

    def add_modules(self, module_list: List[str]):
        """Take a list of modules and create sets of corresponding node/relationship objects

        Args:
            module_list (List[str]): list of socialgene modules (e.g. ["base", "protein"])
        """
        # if only a single module str is provided, make sure it's a lists
        if isinstance(module_list, str):
            module_list = [module_list]
        elif isinstance(module_list, list):
            pass
        else:
            raise ValueError(type(module_list))
        for module in module_list:
            for node in self.modules.get(module).nodes:
                if node_obj := self.nodes.get(node):
                    self.selected_nodes.add(node_obj)
            for rel in self.modules.get(module).relationships:
                if rel_obj := self.relationships.get(rel):
                    self.selected_relationships.add(rel_obj)

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
                outdir, header=node.header, header_filename=node.header_filename
            )
        for rel in self.selected_relationships:
            self._writer(outdir, header=rel.header, header_filename=rel.header_filename)
