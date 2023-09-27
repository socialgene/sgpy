from socialgene.neo4j.schema.relationships import Relationship, Relationships
import pytest


def test_relationship_keys():
    a = Relationships()
    assert list(a.relationships.keys()) == [
        "ANNOTATES",
        "ASSEMBLES_TO",
        "ENCODES",
        "TAXON_PARENT",
        "GO_ANN",
        "PROTEIN_TO_GO",
        "GOTERM_RELS",
        "ROLE_ANN",
        "MAINROLE_ANN",
        "SUBROLE_ANN",
        "IS_TAXON",
        "BLASTP",
        "MMSEQS2",
        "CLUSTER_TO_FILE",
        "MOLECULAR_NETWORK",
        "METABO",
        "SOURCE_DB",
    ]


def test_add_node():
    a = Relationships()
    a.add_relationship(
        neo4j_label="test_neo4j_label",
        description="test_description",
        header_filename="test_header_filename",
        target_subdirectory="test_target_subdirectory",
        target_extension="test_target_extension",
        header=[
            "test_header",
        ],
    )
    assert isinstance(a.relationships["test_neo4j_label"], Relationship)
    assert a.relationships["test_neo4j_label"].__dict__ == {
        "neo4j_label": "test_neo4j_label",
        "description": "test_description",
        "header_filename": "test_header_filename",
        "target_subdirectory": "test_target_subdirectory",
        "target_extension": "test_target_extension",
        "header": ["test_header"],
        "multilabel": False,
    }


def test_add_node_fail():
    a = Relationships()
    with pytest.raises(TypeError):
        a.add_relationship("a")
    with pytest.raises(TypeError):
        a.add_relationship(
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


@pytest.mark.parametrize("rel_label", Relationships().relationships.keys())
def test_node_labels_are_keys(rel_label):
    assert rel_label == Relationships().relationships[rel_label].neo4j_label
