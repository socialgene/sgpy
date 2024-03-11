from socialgene.neo4j.neo4j_element import Node, Relationship
from socialgene.nextflow.nodes import NUCLEOTIDE


class GeneCluster(Node):
    neo4j_label = ["genecluster"]
    description = "A gene cluster"
    property_specifications = {
        "uid": str,
        "start": int,
        "end": int,
    }
    required_properties = ["uid", "start", "end"]


class GeneClusterToNucleotide(Relationship):
    neo4j_label = "ENCODES"
    description = "Connects the nucleotide sequence to the gene cluster"
    property_specifications = {
        "start": int,
        "end": int,
        "core_start": int,
        "core_end": int,
    }
    start_class = NUCLEOTIDE
    end_class = GeneCluster
