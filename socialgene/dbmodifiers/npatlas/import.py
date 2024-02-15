from socialgene.addons.npatlas.nr import NPAtlas, NPAtlasParser, NPAtlasNode, NPAtlasToMibig, NPAtlasToNpmrd,NPAtlasToGnps
import json
import concurrent.futures


def bro(i):
    return NPAtlasParser(i)
    a = NPAtlasParser(i)

def bro2(a):
    return NPAtlasNode(properties={
            "uid": a.uid,
            "original_name": a.original_name,
            "mol_formula": a.mol_formula,
            "mol_weight": a.mol_weight,
            "exact_mass": a.exact_mass,
            "inchikey": a.inchikey,
            "smiles": a.smiles,
            "cluster_id": a.cluster_id,
            "node_id": a.node_id,
            "synonyms": a.synonyms,
            "inchi": a.inchi,
            "m_plus_h": a.m_plus_h,
            "m_plus_na": a.m_plus_na,
            "genus": a.genus,
            "species": a.species,
            "classyfire_class": a.classyfire_class,
            "classyfire_subclass": a.classyfire_subclass,
         })

def process_json(json_path):
    temp1 = NPAtlas(atlas_json_path=json_path)
    list_of_nodes = []
    with open(temp1.path, "r") as f:
        json_data = json.load(f)
        with concurrent.futures.ThreadPoolExecutor() as executor:
            results = executor.map(bro, json_data)
            for result in results:
                list_of_nodes.append(result)
    return list_of_nodes

# Example usage
json_path = "/home/chase/Downloads/ttt/2.json"
list_of_nodes = process_json(json_path)

NPAtlasNode.add_multiple_to_neo4j([bro2(i) for i in list_of_nodes], batch_size=1000, create=True)

np_to_mibig = []
np_to_npmrd = []
np_to_gnps = []
for i in list_of_nodes:
    node=bro2(i)
    if hasattr(i, "mibig"):
        if isinstance(i.mibig, set):
            for x in i.mibig:
                np_to_mibig.append(
                    NPAtlasToMibig(
                        start = x,
                        end=node)
                )
    if hasattr(i, "npmrd"):
        if isinstance(i.npmrd, set):
            for x in i.npmrd:
                np_to_npmrd.append(
                    NPAtlasToNpmrd(
                        start = node,
                        end=x)
                )
    if hasattr(i, "gnps"):
        if isinstance(i.gnps, set):
            for x in i.gnps:
                np_to_gnps.append(
                    NPAtlasToGnps(
                        start = node,
                        end=x)
                )



NPAtlasToMibig.add_multiple_to_neo4j(np_to_mibig, batch_size=1000, create=True)

NPAtlasToNpmrd.add_multiple_to_neo4j(np_to_npmrd, batch_size=1000, create=True)

NPAtlasToGnps.add_multiple_to_neo4j(np_to_gnps, batch_size=1000, create=True)

