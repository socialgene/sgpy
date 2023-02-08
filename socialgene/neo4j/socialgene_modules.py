from socialgene.neo4j.define_nodes import define_nodes
from socialgene.neo4j.define_relationships import define_relationships
from socialgene.neo4j.define_hmm import define_hmms
from socialgene.neo4j.module_abc import Modules

# This module...
# populates the nodes in Modules(); then Modules is accessed as SocialGeneModules()
define_nodes()
define_relationships()
define_hmms()


class SocialGeneModules(Modules):
    pass
