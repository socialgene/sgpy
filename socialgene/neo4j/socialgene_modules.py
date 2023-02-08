from socialgene.neo4j.define_nodes import Nodes
from socialgene.neo4j.define_relationships import Relationships
from socialgene.neo4j.define_hmm import hmm_sources


def parse_hmmlist_input(input):
    # Filter hmm databases based on input list of hmm database names or "all"
    # accept "all" as a list or string
    if input == "all" or "all" in input:
        temp = [i for i in hmm_sources if i != "local"]
        return temp
    else:
        temp = [i for i in input if i in hmm_sources]
        return temp


# auto-generate hmm database nodes and relationship modules
for i in hmm_sources:
    Nodes()._add(
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
    Relationships()._add(
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


class Modules(Nodes, Relationships):
    pass
