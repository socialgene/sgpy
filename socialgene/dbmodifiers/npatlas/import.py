from socialgene.dbmodifiers.classyfire.importer import main as get_classyfire
from socialgene.external_db_classes.classyfire import ClassyFire
from socialgene.external_db_classes.npatlas import Npatlas

get_classyfire()

temp1 = Npatlas.download()

for i in temp1:
    temp = Npatlas(i)
    temp.add_node_to_neo4j()
    temp.merge_with_mibig()
    temp.merge_with_classyfire()

33372
