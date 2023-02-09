from typing import Dict, List, Set
from socialgene.neo4j.schema.define_nodes import Nodes
from socialgene.neo4j.schema.define_relationships import Relationships
from socialgene.neo4j.schema.define_hmmlist import hmm_sources
from socialgene.neo4j.schema.define_modules import Modules
from socialgene.neo4j.schema.node_relationship_class import Neo4jElement


def parse_hmmlist_input(input):
    # Filter hmm databases based on input list of hmm database names or "all"
    # accept "all" as a list or string
    if input == "all" or "all" in input:
        temp = [i for i in hmm_sources if i != "local"]
        return temp
    else:
        temp = [i for i in input if i in hmm_sources]
        return temp


class SocialgeneModules(Relationships, Nodes, Modules):
    def __init__(
        self,
    ):
        super().__init__()
        # auto-generate hmm database nodes and relationship modules
        for i in hmm_sources:
            self.add_node(
                neo4j_label=i,
                header_filename=f"{i}_hmms_out.header",
                target_subdirectory="hmm_tsv_parse",
                target_extension=f"{i}_hmms_out",
                header=[
                    ":IGNORE",
                    "accession",
                    f"id:ID({i})",
                    "description",
                    "category",
                ],
            )
            self.add_relationship(
                neo4j_label="SOURCE_DB",
                header_filename=f"{i}_hmms_out_relationships.header",
                target_subdirectory="hmm_tsv_parse",
                target_extension=f"{i}_hmms_out",
                header=[
                    ":START_ID(hmm)",
                    ":IGNORE",
                    f":END_ID({i})",
                    ":IGNORE",
                    ":IGNORE",
                ],
            )

    def _get_by_label(self, x, y):
        return (i for i in x if i.neo4j_label in y)

    def get_nodes(self, input):
        return self._get_by_label(x=self.nodes, y=input)

    def get_relationships(self, input):
        return self._get_by_label(x=self.relationships, y=input)

    def node_and_rel_dict_by_module_name(
        self, module_list: List[str]
    ) -> Dict[str, Set[Neo4jElement]]:
        """Take a list fo modules and return a dictionary containing a subset of node and relationship Neo4jElements

        Args:
            module_list (List[str]): list of socialgene modules (e.g. "base", "protein", etc)

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
