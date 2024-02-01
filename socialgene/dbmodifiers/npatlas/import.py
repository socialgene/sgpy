from socialgene.addons.npatlas import NPAtlas
from socialgene.dbmodifiers.classyfire.importer import main as get_classyfire

# get_classyfire()

temp1 = NPAtlas("/home/chase/Downloads/NPAtlas_download.json")

temp1._hydrate()


temp1.add_nodes_to_neo4j()
