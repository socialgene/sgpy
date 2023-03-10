from socialgene.neo4j.schema.node_relationship_class import NR


HMM_SOURCES = [
    "amrfinder",
    "antismash",
    "bigslice",
    "classiphage",
    "ipresto",
    "local",
    "pfam",
    "prism",
    "resfams",
    "tigrfam",
    "virus_orthologous_groups",
]


class Hmms(NR):
    def _import(self):
        # auto-generate hmm database nodes and relationship modules
        for i in HMM_SOURCES:
            self.add_node(
                neo4j_label=i,
                header_filename=f"{i}_hmm_source.header",
                target_subdirectory="hmm_tsv_parse",
                target_extension=f"{i}_hmm_source",
                header=[
                    "relative_path",
                    f"name:ID({i})",
                    "acc",
                    "description",
                    "date",
                    ":IGNORE",
                    "model_length",
                    "category",
                    "subcategory",
                ],
            )
            self.add_relationship(
                neo4j_label="SOURCE_DB",
                header_filename=f"{i}_hmm_source_relationships.header",
                target_subdirectory="hmm_tsv_parse",
                target_extension=f"{i}_hmm_source",
                header=[
                    ":IGNORE",
                    f":END_ID({i})",
                    ":IGNORE",
                    ":IGNORE",
                    ":IGNORE",
                    ":START_ID(hmm)",
                    ":IGNORE",
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
            _hmms = HMM_SOURCES
        else:
            _hmms = [i for i in input if i in HMM_SOURCES]
        return _hmms
