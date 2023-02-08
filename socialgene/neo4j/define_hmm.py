from socialgene.neo4j.module_abc import Modules

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


def define_hmms():
    for i in hmm_sources:
        Modules().add_node(
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
        Modules().add_relationship(
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
