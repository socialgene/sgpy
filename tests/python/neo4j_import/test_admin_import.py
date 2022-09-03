from socialgene.neo4j.sg_modules import sgModules, Neo4jImportData


def test_no_overlapping_node_keys():
    a = list(sgModules().nodes.keys())
    assert len(a) == len(set(a))


def test_sgmodules_in_Neo4jImportData():
    # check that all sg_modules actually point to something
    nid = Neo4jImportData()
    sgm = sgModules()
    sgm_list = set([x for xs in sgm.nodes.values() for x in xs])
    sgm_list.update(set([x for xs in sgm.relationships.values() for x in xs]))
    nid_list = set(nid.nodes.keys())
    nid_list.update(set(nid.relationships.keys()))
    assert all(elem in nid_list for elem in sgm_list)
