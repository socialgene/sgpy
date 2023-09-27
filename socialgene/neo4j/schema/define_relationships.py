from socialgene.neo4j.schema.neo4j_element import Neo4jElement
from socialgene.neo4j.schema.node_relationship_class import NR
from rich.console import Console, ConsoleOptions, RenderResult
from rich.table import Table
import builtins
from rich import print
import re

# use rich to print
builtins.print = print


class Relationship(Neo4jElement):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _start(self):
        for i in self.header:
            if "START_ID" in i:
                return i

    def _end(self):
        for i in self.header:
            if "END_ID" in i:
                return i

    @property
    def start(self):
        try:
            return re.findall(r"\((.*?)\)", self._start())[0]
        except Exception as e:
            raise

    @property
    def end(self):
        try:
            return re.findall(r"\((.*?)\)", self._end())[0]
        except Exception as e:
            raise

    @property
    def cypher_string(self):
        return f"(:{self.start})-[:{self.neo4j_label}]->(:{self.end})"

    def __rich_console__(
        self, console: Console, options: ConsoleOptions
    ) -> RenderResult:
        table = Table(title="Relationship")
        table.add_column("Label", justify="left", style="cyan", no_wrap=True, ratio=1)
        table.add_column("Relationship", style="magenta", ratio=1)
        table.add_column("NF results subdirectory", style="magenta", ratio=1)
        table.add_column("Neo4j header file", style="magenta", ratio=1)
        table.add_row(
            self.neo4j_label,
            self.cypher_string,
            self.target_subdirectory,
            self.header_filename,
        )
        yield table


class Relationships:
    __slots__ = ["relationships"]

    def __init__(
        self,
    ):
        self.relationships = set()
        self.add_relationship(
            neo4j_label="ANNOTATES",
            description="",
            header_filename="protein_to_hmm_header.header",
            target_subdirectory="parsed_domtblout",
            target_extension="parseddomtblout",
            header=[
                ":END_ID(protein)",
                ":START_ID(hmm)",
                "env_from:Long",
                "env_to:Long",
                "seq_pro_score:Float",
                "evalue:Long",
                "i_evalue:Long",
                "domain_bias:Float",
                "domain_score:Float",
                "seq_pro_bias:Float",
                "hmm_from:Long",
                "hmm_to:Long",
                "ali_from:Long",
                "ali_to:Long",
                "exponentialized:Boolean",
            ],
        )

        self.add_relationship(
            neo4j_label="ASSEMBLES_TO",
            description="",
            header_filename="assembly_to_locus.header",
            target_subdirectory="genomic_info",
            target_extension="assembly_to_locus",
            header=[":END_ID(assembly)", ":START_ID(nucleotide)"],
        )

        self.add_relationship(
            neo4j_label="ENCODES",
            description="",
            header_filename="locus_to_protein.header",
            target_subdirectory="genomic_info",
            target_extension="locus_to_protein",
            header=[
                ":START_ID(nucleotide)",
                ":END_ID(protein)",
                "protein_id",
                "locus_tag",
                "start:Long",
                "end:Long",
                "strand:Long",
                "description",
                "partial_on_complete_genome:Boolean",
                "missing_start:Boolean",
                "missing_stop:Boolean",
                "internal_stop:Boolean",
                "partial_in_the_middle_of_a_contig:Boolean",
                "missing_N_terminus:Boolean",
                "missing_C_terminus:Boolean",
                "frameshifted:Boolean",
                "too_short_partial_abutting_assembly_gap:Boolean",
                "incomplete:Boolean",
            ],
        )

        self.add_relationship(
            neo4j_label="TAXON_PARENT",
            description="",
            header_filename="taxid_to_taxid.header",
            target_subdirectory="taxdump_process",
            target_extension="taxid_to_taxid",
            header=[":START_ID(taxid)", ":END_ID(taxid)"],
        )

        self.add_relationship(
            neo4j_label="GO_ANN",
            description="",
            header_filename="tigrfam_to_go.header",
            target_subdirectory="tigrfam_info",
            target_extension="tigrfam_to_go",
            header=[":START_ID(hmm_source)", ":END_ID(goterm)"],
        )

        self.add_relationship(
            neo4j_label="PROTEIN_TO_GO",
            description="",
            header_filename="protein_to_go.header",
            target_subdirectory="protein_info",
            target_extension="protein_to_go",
            header=[":START_ID(protein)", ":END_ID(goterm)"],
        )

        self.add_relationship(
            neo4j_label="GOTERM_RELS",
            description="",
            multilabel=True,
            header_filename="go_to_go.header",
            target_subdirectory="goterms",
            target_extension="goterm_edgelist",
            header=[":START_ID(goterm)", ":END_ID(goterm)", ":TYPE"],
        )

        self.add_relationship(
            neo4j_label="ROLE_ANN",
            description="",
            header_filename="tigrfam_to_role.header",
            target_subdirectory="tigrfam_info",
            target_extension="tigrfam_to_role",
            header=[":START_ID(hmm_source)", ":END_ID(tigrfam_role)"],
        )

        self.add_relationship(
            neo4j_label="MAINROLE_ANN",
            description="",
            header_filename="tigrfamrole_to_mainrole.header",
            target_subdirectory="tigrfam_info",
            target_extension="tigrfamrole_to_mainrole",
            header=[":START_ID(tigrfam_role)", ":END_ID(tigrfam_mainrole)"],
        )

        self.add_relationship(
            neo4j_label="SUBROLE_ANN",
            description="",
            header_filename="tigrfamrole_to_subrole.header",
            target_subdirectory="tigrfam_info",
            target_extension="tigrfamrole_to_subrole",
            header=[":START_ID(tigrfam_role)", ":END_ID(tigrfam_subrole)"],
        )

        self.add_relationship(
            neo4j_label="IS_TAXON",
            description="",
            header_filename="assembly_to_taxid.header",
            target_subdirectory="genomic_info",
            target_extension="assembly_to_taxid",
            header=[":START_ID(assembly)", ":END_ID(taxid)"],
        )

        self.add_relationship(
            neo4j_label="BLASTP",
            description="",
            header_filename="blastp.header",
            target_subdirectory="diamond_blastp",
            target_extension="blast6",
            header=[
                ":START_ID(protein)",
                ":END_ID(protein)",
                "pident:Float",
                "length:Long",
                "mismatch:Long",
                "gapopen:Long",
                "qstart:Long",
                "qend:Long",
                "sstart:Long",
                "send:Long",
                "evalue:Float",
                "bitscore:Float",
                "qcovhsp:Float",
            ],
        )

        self.add_relationship(
            neo4j_label="MMSEQS2",
            description="",
            header_filename="mmseqs2.header",
            target_subdirectory="mmseqs2_cluster",
            target_extension="mmseqs2_results_cluster.tsv",
            header=[":START_ID(protein)", ":END_ID(protein)", ":TYPE"],
        )

        self.add_relationship(
            neo4j_label="CLUSTER_TO_FILE",
            description="",
            header_filename="cluster_to_source_file.header",
            target_subdirectory="paired_omics",
            target_extension="cluster_to_source_file",
            header=[":END_ID(mz_cluster_index)", ":START_ID(mz_source_file)"],
        )

        self.add_relationship(
            neo4j_label="MOLECULAR_NETWORK",
            description="",
            header_filename="molecular_network.header",
            target_subdirectory="paired_omics",
            target_extension="molecular_network",
            header=[
                ":START_ID(mz_cluster_index)",
                ":END_ID(mz_cluster_index)",
                "delta_mz:Float",
                "meh:Float",
                "cosine:Float",
                "other_score:Float",
            ],
        )

        self.add_relationship(
            neo4j_label="METABO",
            description="",
            header_filename="assembly_to_mz_file.header",
            target_subdirectory="paired_omics",
            target_extension="assembly_to_mz_file",
            header=[":START_ID(assembly)", ":END_ID(mz_source_file)"],
        )

        self.add_relationship(
            neo4j_label="SOURCE_DB",
            description="",
            header_filename="hmm_source_relationships.header",
            target_subdirectory="hmm_info",
            target_extension=".hmminfo",
            header=[
                ":END_ID(hmm_source)",
                ":IGNORE",
                ":IGNORE",
                ":IGNORE",
                ":IGNORE",
                ":IGNORE",
                ":IGNORE",
                ":IGNORE",
                ":IGNORE",
                ":START_ID(hmm)",
                ":IGNORE",
                ":IGNORE",
                ":IGNORE",
                ":IGNORE",
                ":IGNORE",
                ":IGNORE",
            ],
        )

    def add_relationship(self, **kwargs):
        self.relationships.add(Relationship(**kwargs))

    def __rich_console__(
        self, console: Console, options: ConsoleOptions
    ) -> RenderResult:
        table = Table(title="Relationships", show_lines=True)
        table.add_column("Label", justify="left", style="cyan", no_wrap=True, ratio=1)
        table.add_column("Relationship", style="magenta", ratio=1)
        table.add_column("NF results subdirectory", style="magenta", ratio=1)
        table.add_column("Neo4j header file", style="magenta", ratio=1)
        for i in sorted(list(self.relationships), key=lambda x: x.neo4j_label):
            table.add_row(
                i.neo4j_label,
                i.cypher_string,
                i.target_subdirectory,
                i.header_filename,
            )
        yield table


def print_info():
    print(Relationships())


if __name__ == "__main__":
    print_info()
