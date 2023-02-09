from socialgene.neo4j.socialgene_modules import SocialGeneModules
from socialgene.utils.logging import log
from pathlib import Path
import csv


# `SocialgeneModules` groups `Neo4jImportData`
# `Neo4jImportData` defines information and structure about data that will be imported into Neo4j


hmm_sources = [
    "pfam",
    "antismash",
    "tigrfam",
    "amrfinder",
    "prism",
    "resfams",
    "bigslice",
    "classiphage",
    "virus_orthologous_groups",
    "local",
]


def parse_hmmlist_input(input):
    # Filter hmm databases based on input list of hmm database names or "all"
    # accept "all" as a list or string
    if input == "all" or "all" in input:
        temp = [i for i in hmm_sources if i != "local"]
        return temp
    else:
        temp = [i for i in input if i in hmm_sources]
        return temp


class Neo4jImportData:
    def __init__(
        self,
    ):
        self.nodes = {}

        self.relationships = {}


class SocialgeneModules:
    def __init__(
        self,
    ):
        # if adding both a node and relationship for the same module, give it the same name in both dicts
        self.node_groups = {
            "base": [
                "parameters",
                "assembly",
                "nucleotide",
                "protein",
            ],
            "hmms": [],
            "ncbi_taxonomy": ["taxid"],
            "base_hmm": ["hmm"],
            "tigrfam": [
                "goterm",
                "tigrfam_mainrole",
                "tigrfam_subrole",
                "tigrfam_role",
            ],
            "paired_omics": ["mz_cluster_index", "mz_source_file"],
        }
        self.relationship_groups = {
            "base": [
                "contains",
                "assembles_to",
            ],
            "base_hmm": ["annotates"],
            "hmms": [],
            "tigrfam": [
                "mainrole_ann",
                "role_ann",
                "subrole_ann",
                "go_ann",
            ],
            "ncbi_taxonomy": ["belongs_to", "assembly_to_taxid"],
            "paired_omics": ["cluster_to_file", "molecular_network", "metabo"],
        }

        # enable all node and relationships as individual sgmodules options
        for i in Neo4jImportData().nodes.keys():
            if i not in self.nodes:
                self.nodes[i] = [i]
        for i in Neo4jImportData().relationships.keys():
            if i not in self.relationships:
                self.relationships[i] = [i]
        self.sg_modules = {"nodes": self.nodes, "relationships": self.relationships}

    @staticmethod
    def _filter(input_sg_modules, input_dict):
        temp = {k: v for k, v in input_dict.items() if k in input_sg_modules}
        return temp

    def filter_nodes(self, input_sg_modules):
        if "hmms" in input_sg_modules:
            input_sg_modules.append("base_hmm")
        return self._filter(input_sg_modules, self.nodes)

    def filter_relationships(self, input_sg_modules):
        if "hmms" in input_sg_modules:
            input_sg_modules.append("base_hmm")
        return self._filter(input_sg_modules, self.relationships)


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


def write_neo4j_headers(sg_modules: list, hmmlist: list, outdir: str):
    sg_mod_object = SocialgeneModules()
    header_object = Neo4jImportData()
    reduced_dict = {
        "nodes": sg_mod_object.filter_nodes(sg_modules),
        "relationships": sg_mod_object.filter_relationships(sg_modules),
    }
    for rd_key, rd_value in reduced_dict.items():
        for sg_mod_key, header_keys in rd_value.items():
            if sg_mod_key == "hmms":
                keylist = [i for i in header_keys if i in hmmlist]
                keylist.extend(rd_value["base_hmm"])
            else:
                keylist = header_keys
            for i in keylist:
                _writer(outdir, getattr(header_object, rd_key)[i])
