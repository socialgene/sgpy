from dataclasses import dataclass
from typing import List

import socialgene.nextflow.nodes as nodes
import socialgene.nextflow.relationships as rels


@dataclass
class Single_Module:
    """Class for keeping track of an item in inventory."""

    module_id: str
    nodes: List[str]
    relationships: List[str]


class Modules:
    # These are used in Nextflow, here: https://github.com/socialgene/sgnf/blob/main/subworkflows/local/sg_modules.nf
    # They simply group the different nodes/relationships into somewhat sensible "modules"
    def __init__(self, *args, **kwargs):
        self.modules = {}
        self._add_module(
            module_id="base",
            nodes=[nodes.ASSEMBLY, nodes.NUCLEOTIDE, nodes.PROTEIN],
            relationships=[
                rels.ASSEMBLES_TO,
                rels.ENCODES,
            ],
        )
        self._add_module(
            module_id="go",
            nodes=[nodes.GOTERM],
            relationships=[
                rels.PROTEIN_TO_GO,
                rels.GOTERM_RELS,
            ],
        )
        self._add_module(
            module_id="protein",
            nodes=[nodes.PROTEIN],
            relationships=[],
        )
        self._add_module(
            module_id="parameters",
            nodes=[nodes.PARAMETERS],
            relationships=[],
        )
        self._add_module(
            module_id="base_hmm",
            nodes=[
                nodes.HMM,
                nodes.HMM_SOURCE,
            ],
            relationships=[rels.ANNOTATES, rels.SOURCE_DB],
        )
        self._add_module(
            module_id="ncbi_taxonomy",
            nodes=[nodes.TAXID],
            relationships=[rels.TAXON_PARENT, rels.IS_TAXON],
        )
        self._add_module(
            module_id="tigrfam",
            nodes=[
                nodes.GOTERM,
                nodes.TIGRFAM_MAINROLE,
                nodes.TIGRFAM_SUBROLE,
                nodes.TIGRFAM_ROLE,
            ],
            relationships=[
                rels.MAINROLE_ANN,
                rels.ROLE_ANN,
                rels.SUBROLE_ANN,
                rels.GO_ANN,
            ],
        )
        self._add_module(
            module_id="blastp",
            nodes=[],
            relationships=[rels.BLASTP],
        )
        self._add_module(
            module_id="mmseqs",
            nodes=[],
            relationships=[rels.MMSEQS2],
        )

    def _add_module(self, **kwargs):
        _module_id = kwargs.get("module_id")
        if _module_id in self.modules:
            raise ValueError(f"Module ID: {_module_id} already exists")
        self.modules.update({_module_id: Single_Module(**kwargs)})
