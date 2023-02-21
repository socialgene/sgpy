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
    def _add_module(self, **kwargs):
        _module_id = kwargs.get("module_id")
        self.modules.update({_module_id: Single_Module(**kwargs)})

    def __init__(
        self,
    ):
        self.modules = {}
        self._add_module(
            module_id="base",
            nodes=[
                "assembly",
                "nucleotide",
                "protein",
            ],
            relationships=[
                "CONTAINS",
                "ASSEMBLES_TO",
                "ENCODES",
            ],
        )
        self._add_module(
            module_id="protein",
            nodes=["protein"],
            relationships=[],
        )
        self._add_module(
            module_id="parameters",
            nodes=["parameters"],
            relationships=[],
        )
        self._add_module(
            module_id="base_hmm",
            nodes=["hmm"],
            relationships=["ANNOTATES"],
        )
        self._add_module(
            module_id="ncbi_taxonomy",
            nodes=["taxid"],
            relationships=["BELONGS_TO", "ASSEMBLY_TO_TAXID"],
        )
        self._add_module(
            module_id="tigrfam",
            nodes=[
                "goterm",
                "tigrfam_mainrole",
                "tigrfam_subrole",
                "tigrfam_role",
            ],
            relationships=[
                "MAINROLE_ANN",
                "ROLE_ANN",
                "SUBROLE_ANN",
                "GO_ANN",
            ],
        )
        self._add_module(
            module_id="blastp",
            nodes=[],
            relationships=["BLASTP"],
        )
        self._add_module(
            module_id="mmseqs",
            nodes=[],
            relationships=["MMSEQS2"],
        )
        self._add_module(
            module_id="paired_omics",
            nodes=["mz_cluster_index", "mz_source_file"],
            relationships=["CLUSTER_TO_FILE", "MOLECULAR_NETWORK", "METABO"],
        )
