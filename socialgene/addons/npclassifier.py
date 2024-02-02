from socialgene.neo4j.neo4j_element import Node, Relationship


class NPClassifierClass(Node):
    def __init__(self, *args, **kwargs):
        super().__init__(
            neo4j_label="npclassifier_class",
            description="Represents a NPClassifier class",
            properties={
                "uid": "string",
            },
        )


class NPClassifierPathway(Node):
    def __init__(self, *args, **kwargs):
        super().__init__(
            neo4j_label="npclassifier_pathway",
            description="Represents a NPClassifier pathway",
            properties={
                "uid": "string",
            },
        )


class NPClassifierSuperclass(Node):
    def __init__(self, *args, **kwargs):
        super().__init__(
            neo4j_label="npclassifier_superclass",
            description="Represents a NPClassifier superclass",
            properties={
                "uid": "string",
            },
        )


class NPClassifierIsA(Relationship):
    def __init__(self, *args, **kwargs):
        super().__init__(
            neo4j_label="IS_A",
            description="Represents a relationship between npclassifier nodes",
            start=NPClassifierClass,
            end=NPClassifierClass,
        )
