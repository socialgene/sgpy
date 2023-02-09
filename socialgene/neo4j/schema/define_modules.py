from typing import List
from dataclasses import dataclass


@dataclass
class Single_Module:
    """Class for keeping track of an item in inventory."""

    module_id: str
    nodes: List[str]
    relationships: List[str]


class Modules:
    # These are used in Nextflow, here: https://github.com/socialgene/sgnf/blob/main/subworkflows/local/sg_modules.nf
    # They simply group the different files/inputs into hopefully-sensible "modules"
    # It could be coded differently, but this way will hopefully make it easier for someone to add on a group
    def add_module(self, **kwargs):
        _module_id = kwargs.get("module_id")
        self.modules.update({_module_id: Single_Module(**kwargs)})

    def __init__(
        self,
    ):
        self.modules = {}
        self.add_module(
            module_id="base",
            nodes=[
                "parameters",
                "assembly",
                "nucleotide",
                "protein",
            ],
            relationships=[
                "contains",
                "assembles_to",
            ],
        )
        self.add_module(
            module_id="base_hmm",
            nodes=["hmm"],
            relationships=["annotates"],
        )
        self.add_module(
            module_id="ncbi_taxonomy",
            nodes=["taxid"],
            relationships=["belongs_to", "assembly_to_taxid"],
        )
        self.add_module(
            module_id="tigrfam",
            nodes=[
                "goterm",
                "tigrfam_mainrole",
                "tigrfam_subrole",
                "tigrfam_role",
            ],
            relationships=[
                "mainrole_ann",
                "role_ann",
                "subrole_ann",
                "go_ann",
            ],
        )
        self.add_module(
            module_id="paired_omics",
            nodes=["mz_cluster_index", "mz_source_file"],
            relationships=["cluster_to_file", "molecular_network", "metabo"],
        )
