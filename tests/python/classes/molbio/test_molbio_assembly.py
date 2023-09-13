from collections import OrderedDict
from unittest.mock import patch

import pytest

from socialgene.base.molbio import Assembly, Taxonomy


def test_fail1():
    with pytest.raises(TypeError):
        Assembly()


def test_id():
    a = Assembly(id="hi")
    assert a.id == "hi"


def test_name():
    a = Assembly(id="hi")
    assert a.name == "hi"


def test_dict():
    a = Assembly(id="hi")
    b = a.__dict__
    assert isinstance(b["taxonomy"], Taxonomy)
    del b["taxonomy"]
    assert b == {
        "id": "hi",
        "info": OrderedDict(
            [
                ("mol_type", None),
                ("altitude", None),
                ("bio_material", None),
                ("bioproject", None),
                ("biosample", None),
                ("cell_line", None),
                ("cell_type", None),
                ("chromosome", None),
                ("clone", None),
                ("clone_lib", None),
                ("collected_by", None),
                ("collection_date", None),
                ("country", None),
                ("cultivar", None),
                ("culture_collection", None),
                ("db_xref", None),
                ("dev_stage", None),
                ("ecotype", None),
                ("environmental_sample", None),
                ("focus", None),
                ("germline", None),
                ("haplogroup", None),
                ("haplotype", None),
                ("host", None),
                ("identified_by", None),
                ("isolate", None),
                ("isolation_source", None),
                ("lab_host", None),
                ("lat_lon", None),
                ("macronuclear", None),
                ("map", None),
                ("mating_type", None),
                ("metagenome_source", None),
                ("note", None),
                ("organelle", None),
                ("PCR_primers", None),
                ("plasmid", None),
                ("pop_variant", None),
                ("proviral", None),
                ("rearranged", None),
                ("segment", None),
                ("serotype", None),
                ("serovar", None),
                ("sex", None),
                ("specimen_voucher", None),
                ("strain", None),
                ("sub_clone", None),
                ("submitter_seqid", None),
                ("sub_species", None),
                ("sub_strain", None),
                ("tissue_lib", None),
                ("tissue_type", None),
                ("transgenic", None),
                ("type_material", None),
                ("variety", None),
            ]
        ),
        "loci": {},
        "name": "hi",
        "taxid": None,
    }


@patch("socialgene.base.molbio.uuid4")
def test_add_locus_uuid(my_method):
    my_method.return_value = "ASDSDSD"
    a = Assembly(id="a")
    a.add_locus()
    assert list(a.loci.keys())[0] == "ASDSDSD"


def test():
    a = Assembly(id="a")
    a.add_locus("b")
    a.loci["b"].add_feature(start=1000, end=2000)
    a.loci["b"].add_feature(start=1000, end=2000, protein_hash="hiya")
    assert len(a.protein_hash_set) == 2
    assert None in a.protein_hash_set
    assert "hiya" in a.protein_hash_set
    z = a.get_all_proteins()
    assert len(z) == 2
    assert None in z
    assert "hiya" in z


def test_sort_by_middle():
    a = Assembly(id="a")
    a.add_locus("b")
    a.loci["b"].add_feature(start=100, end=5000)
    a.loci["b"].add_feature(start=1000, end=2000, protein_hash="hiya")
    a = [i for i in a.loci["b"].features_sorted_by_midpoint]
    assert a[0].end == 2000
    assert a[1].end == 5000
