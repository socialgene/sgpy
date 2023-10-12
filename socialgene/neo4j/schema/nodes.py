import builtins

from rich import print
from rich.console import Console, ConsoleOptions, RenderResult
from rich.table import Table

from socialgene.base.molbio import LocusAssemblyMetadata
from socialgene.neo4j.schema.neo4j_element import Neo4jElement
from socialgene.utils.lists_to_markdown import markdown_table_from_list

# use rich to print
builtins.print = print


class Node(Neo4jElement):
    """Represents a single Node"""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __rich_console__(
        self, console: Console, options: ConsoleOptions
    ) -> RenderResult:  # pragma: no cover
        table = Table(title="Node")
        table.add_column("Label", justify="left", style="cyan", no_wrap=True, ratio=1)
        table.add_column(
            "Description",
            justify="left",
            style="cyan",
            no_wrap=False,
            ratio=4,
            max_width=50,
        )
        table.add_column("Nextflow results subdirectory", style="magenta", ratio=1)
        table.add_column("Neo4j header file", style="magenta", ratio=1)
        table.add_row(
            self.neo4j_label,
            self.target_subdirectory,
            self.target_subdirectory,
            self.header_filename,
        )
        yield table


class Nodes:
    """Represents multiple Nodes and is where all non-addon Nodes are defined"""

    def __init__(
        self,
    ):
        super().__init__()
        self.nodes = {}
        self.add_node(
            neo4j_label="parameters",
            description="Parameters and environmental variables used during database creation",
            header_filename="parameters.header",
            target_subdirectory="parameters",
            target_extension="socialgene_parameters",
            header=[
                "uid:ID(when)",
                "SG_LOC_NEO4J",
                "SG_LOC_HMMS",
                "NEO4J_dbms_memory_pagecache_size",
                "NEO4J_dbms_memory_heap_initial__size",
                "NEO4J_dbms_memory_heap_max__size",
                "HMMSEARCH_IEVALUE",
                "HMMSEARCH_BACKGROUND",
                "HMMSEARCH_BIASFILTER",
                "HMMSEARCH_NULL2",
                "HMMSEARCH_SEED",
                "HMMSEARCH_Z",
                "HMMSEARCH_DOMZ",
                "HMMSEARCH_F1",
                "HMMSEARCH_F2",
                "HMMSEARCH_F3",
                "HMMSEARCH_E",
                "HMMSEARCH_DOME",
                "HMMSEARCH_INCE",
                "HMMSEARCH_INCDOME",
                "HMMSEARCH_BITCUTOFFS",
                "platform",
                "architecture",
                "py_executable",
                "py_version",
                "genome_download_command",
            ],
        )

        self.add_node(
            neo4j_label="assembly",
            description="Represents a single genome/assembly/BGC. If the input was a FASTA file or if assembly wasn't in the genbank metadata then this will represent the file the data came from.",
            header_filename="assembly.header",
            target_subdirectory="genomic_info",
            target_extension="assemblies",
            header=["uid:ID(assembly)"] + sorted(LocusAssemblyMetadata.__slots__),
        )

        self.add_node(
            neo4j_label="nucleotide",
            description="Represents a single nucleotide sequence (e.g. a contig/scaffold/chromosome)",
            header_filename="locus.header",
            target_subdirectory="genomic_info",
            target_extension="loci",
            header=["uid:ID(nucleotide)"]
            + ["external_id"]
            + LocusAssemblyMetadata.__slots__,
        )

        self.add_node(
            neo4j_label="protein",
            description="Represents a non-redundant protein",
            header_filename="protein_ids.header",
            target_subdirectory="protein_info",
            target_extension="protein_ids",
            header=[
                "uid:ID(protein)",
                "crc64",
            ],
        )
        self.add_node(
            neo4j_label="goterm",
            description="Represent a GO term",
            header_filename="goterms.header",
            target_subdirectory="goterms",
            target_extension="goterms",
            header=["uid:ID(goterm)", "name", "namespace", "def"],
        )

        self.add_node(
            neo4j_label="tigrfam_role",
            description="Represents a TIGRFAM role",
            header_filename="tigrfam_role.header",
            target_subdirectory="tigrfam_info",
            target_extension="tigrfam_role",
            header=["uid:ID(tigrfam_role)"],
        )
        self.add_node(
            neo4j_label="tigrfam_mainrole",
            description="Represents a TIGRFAM main role",
            header_filename="tigrfam_mainrole.header",
            target_subdirectory="tigrfam_info",
            target_extension="tigrfam_mainrole",
            header=["uid:ID(tigrfam_mainrole)"],
        )

        self.add_node(
            neo4j_label="tigrfam_subrole",
            description="Represents a TIGRFAM sub role",
            header_filename="tigrfam_subrole.header",
            target_subdirectory="tigrfam_info",
            target_extension="tigrfam_subrole",
            header=["uid:ID(tigrfam_subrole)"],
        )

        self.add_node(
            neo4j_label="taxid",
            description="Represents a single taxon within NCBI taxonomy",
            header_filename="taxid.header",
            target_subdirectory="taxdump_process",
            target_extension="nodes_taxid",
            header=["uid:ID(taxid)", "name", "rank"],
        )

        self.add_node(
            neo4j_label="mz_cluster_index",
            description="Represents a single m/z cluster",
            header_filename="mz_cluster_index_nodes.header",
            target_subdirectory="paired_omics",
            target_extension="mz_cluster_index_nodes",
            header=[
                "uid:ID(mz_cluster_index)",
                "component_index:Long",
                "parent_mass:Float",
                "precursor_mass:Float",
                "sum_precursor_intensity:Float",
                "Smiles:String",
                "rt_mean::Float",
                "rt_std_err:Float",
                "library_id:String",
                "mq_score:Float",
                "mz_error_ppm:Float",
                "mass_diff:String",
            ],
        )

        self.add_node(
            neo4j_label="mz_source_file",
            description="Represents the file m/z features came from",
            header_filename="mz_source_file.header",
            target_subdirectory="paired_omics",
            target_extension="hash.mz_source_file",
            header=[
                "uid:ID(mz_source_file)",
            ],
        )

        self.add_node(
            neo4j_label="hmm_source",
            description="Represents the source of an HMM model (e.g. PFAM)",
            header_filename="hmm_source.header",
            target_subdirectory="hmm_info",
            target_extension="hmminfo",
            header=[
                "uid:ID(hmm_source)",
                ":LABEL",
                "rel_path:String",
                "name:String",
                "acc:String",
                "notes:String",
                "description:String",
                "date:String",
                "hash:String",
                "hash_used:String",
                "model_length:String",
                "category:String",
                "subcategory:String",
                "ga:String",
                "tc:String",
                "nc:String",
            ],
        )

        self.add_node(
            neo4j_label="hmm",
            description="Represents a single non-redundant HMM model",
            header_filename="sg_hmm_nodes.header",
            target_subdirectory="hmm_info",
            target_extension="sg_hmm_nodes",
            header=["uid:ID(hmm)"],
        )

    def add_node(self, neo4j_label, **kwargs):
        self.nodes[neo4j_label] = Node(neo4j_label=neo4j_label, **kwargs)

    def __rich_console__(
        self, console: Console, options: ConsoleOptions
    ) -> RenderResult:  # pragma: no cover
        table = Table(title="Nodes", show_lines=True)
        table.add_column("Label", justify="left", style="cyan", no_wrap=True, ratio=1)
        table.add_column(
            "Description",
            justify="left",
            style="cyan",
            no_wrap=False,
            ratio=4,
            max_width=50,
        )
        table.add_column("NF results subdirectory", style="magenta", ratio=1)
        table.add_column("Neo4j header file", style="magenta", ratio=1)
        # sort by label which is the key
        for i in (self.nodes[i] for i in sorted(self.nodes.keys())):
            table.add_row(
                i.neo4j_label, i.description, i.target_subdirectory, i.header_filename
            )
        yield table

    def _markdown_table(self):
        cols = [
            (
                "Label",
                "Description",
                "NF results subdirectory",
                "Neo4j header file",
            )
        ]
        rows = [
            (
                i.neo4j_label,
                i.description,
                i.target_subdirectory,
                i.header_filename,
            )
            for i in (self.nodes[i] for i in sorted(self.nodes.keys()))
        ]
        cols.extend(rows)
        print(
            markdown_table_from_list(
                cols,
                align="left",
            )
        )


def print_info():  # pragma: no cover
    print(Nodes())


def print_markdown():  # pragma: no cover
    Nodes()._markdown_table()


def printer():  # pragma: no cover
    import argparse

    parser = argparse.ArgumentParser(description="Print node info")
    parser.add_argument(
        "--markdown",
        help="",
        default=False,
        required=False,
        action=argparse.BooleanOptionalAction,
    )
    args = parser.parse_args()
    if args.markdown:
        print_markdown()
    else:
        print_info()


if __name__ == "__main__":
    print_info()
