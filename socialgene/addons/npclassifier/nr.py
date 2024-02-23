from socialgene.neo4j.neo4j_element import Node


class NPClassifierClass(Node):
    neo4j_label = ["npclassifier_class"]
    description = "Represents a NPClassifier class"
    property_specification = {
        "uid": str,
    }
    constraints_unique = ["uid"]


class NPClassifierPathway(Node):
    neo4j_label = ["npclassifier_pathway"]
    description = "Represents a NPClassifier pathway"
    property_specification = {
        "uid": str,
    }
    constraints_unique = ["uid"]


class NPClassifierSuperclass(Node):
    neo4j_label = ["npclassifier_superclass"]
    description = "Represents a NPClassifier superclass"
    property_specification = {
        "uid": str,
    }
    constraints_unique = ["uid"]
