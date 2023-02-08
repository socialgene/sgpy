from socialgene.neo4j.define_nodes import define_nodes
from socialgene.neo4j.define_relationships import define_relationships
from socialgene.neo4j.module_abc import Modules
from rich import inspect


define_nodes()
define_relationships()
inspect(Modules)
print(len(Modules.nodes))
print(len(Modules.relationships))
