from socialgene.neo4j.schema.node_relationship_class import NR


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


class Hmms(NR):
    def _import(self):
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

    def get_relationships(self, source):
        return (i for i in self.relationships if i.target_extension.startswith(source))

    @staticmethod
    def parse_hmmlist_input(input):
        # Filter hmm databases based on input list of hmm database names or "all"
        # accept "all" as a list or string
        if isinstance(input, str):
            input = [input]
        if "all" in input:
            _hmms = hmm_sources
        else:
            _hmms = [i for i in input if i in hmm_sources]
        return _hmms
