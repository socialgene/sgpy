from socialgene.base.socialgene import SocialGene
from socialgene.clustermap.clustermap import Clustermap, ClustermapUuids, UuidCount


def test_UuidCount_slots():
    a = UuidCount()
    assert a.__slots__ == ["uuid_counter"]


def test_UuidCount_initial_value():
    a = UuidCount()
    assert a.uuid_counter == 0
    b = UuidCount()
    a.uuid_counter += 1
    assert a.uuid_counter == 1
    # b should still be zero
    assert b.uuid_counter == 0
    b.uuid_counter += 10
    assert b.uuid_counter == 10
    assert a.uuid_counter == 1


def test_ClustermapUuids_slots():
    a = ClustermapUuids()
    assert a.__slots__ == ["uuid_counter"]


def test_ClustermapUuids_initial_value():
    a = ClustermapUuids()
    assert a.uuid_counter == 0
    b = ClustermapUuids()
    a.uuid_counter += 1
    assert a.uuid_counter == 1
    # b should still be zero
    assert b.uuid_counter == 0
    b.uuid_counter += 10
    assert b.uuid_counter == 10
    assert a.uuid_counter == 1


def test_build_protein_key():
    temp = ClustermapUuids().build_protein_key(
        assembly_key="assembly-key", locus_key="locus-key", locus_value="locus-value"
    )
    assert temp == "assembly-key_locus-key_locus-value"


def test_empty_sg():
    cm_object = Clustermap()
    sg_object = SocialGene()
    cm_object.create_clustermap_uuids(sg_object=sg_object)
    cm_object.add_cluster(sg_object=sg_object)
    cm_object.add_groups(sg_object=sg_object)
    cm_object.add_links(sg_object=sg_object)
    assert cm_object.primary_assembly is None
    assert cm_object.clusters == []
    assert cm_object.groups == []
    assert cm_object.links == []
    assert cm_object.prot_loc == {}


def test_sg():
    cm_object = Clustermap()
    sg_object = SocialGene()
    _ = sg_object.add_protein(hash_id="assembly_1_locus_1_protein_1")
    sg_object.add_assembly("assembly_1")
    sg_object.assemblies["assembly_1"].add_locus("assembly_1_locus_1")
    sg_object.assemblies["assembly_1"].loci["assembly_1_locus_1"]
    sg_object.assemblies["assembly_1"].loci["assembly_1_locus_1"].add_feature(
        type="protein",
        protein_hash="assembly_1_locus_1_protein_1",
        start=1,
        end=10,
        strand=1,
    )
    cm_object.create_clustermap_uuids(sg_object=sg_object)
    cm_object.add_cluster(sg_object=sg_object)
    cm_object.add_groups(sg_object=sg_object)
    cm_object.add_links(sg_object=sg_object)
    assert cm_object.primary_assembly is None
    assert cm_object.clusters == [
        {
            "uid": "uuid_0",
            "name": "assembly_1",
            "loci": [
                {
                    "uid": "uuid_1",
                    "name": "assembly_1_locus_1",
                    "start": 1,
                    "end": 10,
                    "genes": [
                        {
                            "uid": "uuid_2",
                            "label": None,
                            "names": {"name": None, "description": None},
                            "start": 1,
                            "end": 10,
                            "strand": 1,
                        }
                    ],
                }
            ],
        }
    ]
    assert cm_object.groups == []
    assert cm_object.links == []
    assert cm_object.prot_loc == {
        "assembly_1_locus_1_protein_1": [
            "assembly_1_assembly_1_locus_1_assembly_1_locus_1_protein_1"
        ]
    }
