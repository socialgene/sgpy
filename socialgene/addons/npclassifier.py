from socialgene.neo4j.neo4j_element import Node, Relationship


class NPClassifierClass(Node):
    neo4j_label = "npclassifier_class"
    description = "Represents a NPClassifier class"
    property_specification = {
        "uid": str,
    }


class NPClassifierPathway(Node):
    neo4j_label = "npclassifier_pathway"
    description = "Represents a NPClassifier pathway"
    property_specification = {
        "uid": str,
    }


class NPClassifierSuperclass(Node):
    neo4j_label = "npclassifier_superclass"
    description = "Represents a NPClassifier superclass"
    property_specification = {
        "uid": str,
    }


class NPClassifierIsA(Relationship):
    neo4j_label = "IS_A"
    description = "Represents a relationship between npclassifier nodes"
    start_class = NPClassifierClass
    end_class = NPClassifierClass
