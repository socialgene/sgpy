import collections

from socialgene.neo4j.schema.graph_schema import GraphSchema

a = GraphSchema()


def test_for_node_duplicates():
    """Test for duplicate Neo4j node labels in the schema."""
    node_labels = [i.neo4j_label for i in a.ALL_NODES]
    duplicated = [
        item for item, count in collections.Counter(node_labels).items() if count > 1
    ]
    def_dict = collections.defaultdict(list)
    for i in a.ALL_NODES:
        if i.neo4j_label in duplicated:
            def_dict[i.neo4j_label].append(i.__module__)

    if len(node_labels) != len(set(node_labels)):
        raise ValueError(f"Duplicate node label definitions: {dict(def_dict)}")


# relationship label duplicates are fine
# def test_for_relationship_duplicates():
#     relationship_labels = [i()._Neo4jElement__neo4j_label for i in a.ALL_RELATIONSHIPS]
#     duplicated = [
#         item
#         for item, count in collections.Counter(relationship_labels).items()
#         if count > 1
#     ]
#     def_dict = collections.defaultdict(list)
#     for i in a.ALL_RELATIONSHIPS:
#         if i()._Neo4jElement__neo4j_label in duplicated:
#             def_dict[i()._Neo4jElement__neo4j_label].append(i.__module__)

#     if len(relationship_labels) != len(set(relationship_labels)):
#         raise ValueError(f"Duplicate relationship label definitions: {dict(def_dict)}")
