from socialgene.neo4j.schema.nodes import Node, Nodes
import pytest


def test_node_keys():
    a = Nodes()
    assert list(a.nodes.keys()) == [
        "parameters",
        "assembly",
        "nucleotide",
        "protein",
        "goterm",
        "tigrfam_role",
        "tigrfam_mainrole",
        "tigrfam_subrole",
        "taxid",
        "mz_cluster_index",
        "mz_source_file",
        "hmm_source",
        "hmm",
    ]


def test_add_node():
    a = Nodes()
    a.add_node(
        neo4j_label="test_neo4j_label",
        description="test_description",
        header_filename="test_header_filename",
        target_subdirectory="test_target_subdirectory",
        target_extension="test_target_extension",
        header=[
            "test_header",
        ],
    )
    assert isinstance(a.nodes["test_neo4j_label"], Node)
    assert a.nodes["test_neo4j_label"].__dict__ == {
        "neo4j_label": "test_neo4j_label",
        "description": "test_description",
        "header_filename": "test_header_filename",
        "target_subdirectory": "test_target_subdirectory",
        "target_extension": "test_target_extension",
        "header": ["test_header"],
        "multilabel": False,
    }


def test_add_node_fail():
    a = Nodes()
    with pytest.raises(TypeError):
        a.add_node("a")
    with pytest.raises(TypeError):
        a.add_node(
            neo4j_label="test_neo4j_label",
            description="test_description",
            header_filename="test_header_filename",
            target_subdirectory="test_target_subdirectory",
            target_extension="test_target_extension",
            header=[
                "test_header",
            ],
            bad="a",
        )


@pytest.mark.parametrize("node_label", Nodes().nodes.keys())
def test_node_labels_are_keys(node_label):
    assert node_label == Nodes().nodes[node_label].neo4j_label
