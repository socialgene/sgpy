from socialgene.addons.npatlas import Npatlas
from socialgene.dbmodifiers.classyfire.importer import main as get_classyfire

get_classyfire()

temp1 = Npatlas.download()

for i in temp1:
    temp = Npatlas(i)
    temp.add_node_to_neo4j()
    temp.merge_with_mibig()
    temp.merge_with_classyfire()
