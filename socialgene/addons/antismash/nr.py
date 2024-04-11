from socialgene.neo4j.neo4j_element import Node, Relationship
from socialgene.nextflow.nodes import NUCLEOTIDE





class Product(Node):
    neo4j_label = ["product"]
    description = "Represents the product of an antiSMASH BGC"
    property_specification = {
        "uid": int,
    }
    constraints_unique = ["uid"]


class Category(Node):
    neo4j_label = ["category"]
    description = "Represents the category of an antiSMASH BGC"
    property_specification = {
        "uid": int,
    }
    constraints_unique = ["uid"]

class ProductToCategory(Relationship):
    neo4j_label = "IS_A"
    description = "Connects an antiSMASH product to category "
    start_class = Product
    end_class = Category


class GeneClusterToProduct(Relationship):
    neo4j_label = "IS_A"
    description = "Connects an antiSMASH BGC to product "
    start_class = Product
    end_class = Category
    property_specification = {
        "tool": str,
        "start":int,
        "end":int,
        "core_start":int,
        "core_end":int,
    }

