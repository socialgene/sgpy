from socialgene.neo4j.neo4j_element import Node


class SPECTRUM(Node):
    def __init__(self, *args, **kwargs):
        super().__init__(
            neo4j_label="spectrum",
            description="Represents a GNPS molecular networking spectrum",
            properties={
                "uid": "string",
                "original_filename": "string",
                "parentmass": "float",
                "charge": "int",
                "rettime": "float",
                "assembly": "string",
            },
        )


class NPCLASSIFIER_CLASS(Node):
    def __init__(self, *args, **kwargs):
        super().__init__(
            neo4j_label="npclassifier_class",
            description="Represents a NPClassifier class",
            properties={
                "uid": "string",
            },
        )


class NPCLASSIFIER_PATHWAY(Node):
    def __init__(self, *args, **kwargs):
        super().__init__(
            neo4j_label="npclassifier_pathway",
            description="Represents a NPClassifier pathway",
            properties={
                "uid": "string",
            },
        )


class NPCLASSIFIER_SUPERCLASS(Node):
    def __init__(self, *args, **kwargs):
        super().__init__(
            neo4j_label="npclassifier_superclass",
            description="Represents a NPClassifier superclass",
            properties={
                "uid": "string",
            },
        )
