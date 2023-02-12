import csv
from pathlib import Path
from typing import Dict, List, Set
from socialgene.neo4j.schema.define_nodes import Nodes
from socialgene.neo4j.schema.define_relationships import Relationships
from socialgene.neo4j.schema.define_hmmlist import Hmms, hmm_sources
from socialgene.neo4j.schema.define_modules import Modules
from socialgene.neo4j.schema.node_relationship_class import Neo4jElement
from socialgene.utils.logging import log


class SocialgeneModules(Modules):
    def __init__(
        self,
    ):
        super().__init__()
        self.nodes = set()
        self.relationships = set()
        self.all_hmms = Hmms()
        self.all_nodes = Nodes()
        self.all_relationships = Relationships()

    def add_hmms(self, hmm_list):
        _hmms = []
        if isinstance(hmm_list, str):
            hmm_list = [hmm_list]
        #  if 'all', use all hmms, otherwise filter based on input list
        if "all" in hmm_list:
            _hmms = hmm_sources
        else:
            _hmms = [i for i in hmm_list if i in hmm_sources]
        # add nodes and relationships
        if _hmms:
            for node in self.all_hmms.get_nodes(_hmms):
                self.nodes.add(node)
            for source in _hmms:
                for rel in self.all_hmms.get_relationships(source):
                    self.relationships.add(rel)

    def add_modules(self, module_list: List[str]):
        if isinstance(module_list, str):
            module_list = [module_list]

        for module in module_list:
            module_def = self.modules.get(module)
            for node in self.all_nodes.get_nodes(module_def.nodes):
                self.nodes.add(node)
            for rel in list(
                self.all_relationships.get_relationships(module_def.relationships)
            ):
                self.relationships.add(rel)

    def make_node_and_rel_dict_by_module_name(
        self,
        module_list: List[str],
        hmm_list: List[str],
    ) -> Dict[str, Set[Neo4jElement]]:
        """Take a list fo modules and return a dictionary containing a subset of node and relationship Neo4jElements

        Args:
            module_list (List[str]): list of socialgene modules (e.g. ["base", "protein"])
            hmm_list (List[str]): list of HMM sources (e.g. ["antismash", "pfam"])

        Returns:
            Dict[str, Set[Neo4jElement]]: dictionary containing a subset of node and relationship Neo4jElements
        """
        to_return = {"nodes": set(), "relationships": set()}

        for module in module_list:
            module_def = self.modules.get(module)
            for n in self.get_nodes(module_def.nodes):
                to_return["nodes"].add(n)
            for r in self.get_nodes(module_def.nodes):
                to_return["relationships"].add(r)
        return to_return

    def _writer(outdir, header_dict):
        outpath = Path(outdir, f"{header_dict['header_filename']}")
        with open(outpath, "w") as tsv_output_con:
            tsv_writer = csv.writer(
                tsv_output_con,
                delimiter="\t",
                quotechar='"',
                quoting=csv.QUOTE_MINIMAL,
            )
            tsv_writer.writerow(header_dict["header"])
            log.info(f"\tWriting {header_dict['header_filename']} to: {outpath}")


def write_neo4j_headers(self, module_list: list, hmm_list: list, outdir: str):
    reduced_dict = self.make_node_and_rel_dict_by_module_name(
        module_list=module_list, hmm_list=hmm_list
    )
    for rd_key, rd_value in reduced_dict.items():
        for sg_mod_key, header_keys in rd_value.items():
            for i in keylist:
                _writer(outdir, getattr(header_object, rd_key)[i])
