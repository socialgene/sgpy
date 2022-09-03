import os
import hashlib
from socialgene.parsers.hmm import HMMParser

FIXTURE_DIR = os.path.dirname(os.path.realpath(__file__))
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)
FIXTURE_DIR = os.path.join(FIXTURE_DIR, "data")

pks = os.path.join(FIXTURE_DIR, "pks.hmm")


def test_hmms_object_structure():
    hmms_object = HMMParser()
    assert hmms_object.category is None
    assert hmms_object.hmm_dbs == [
        "prism",
        "bigslice",
        "antismash",
        "amrfinder",
        "resfams",
        "tigrfam",
        "pfam",
        "classiphage",
        "virus_orthologous_groups",
        "local",
    ]
    assert hmms_object.model_info_dict == {}
    assert hmms_object.model_source is None
    assert hmms_object.model_text_dict == {}
    assert hmms_object.nr_models == {}
    assert hmms_object.rel_path is None
    assert hmms_object.single_model_dict == {
        "source": None,
        "source_id": None,
        "rel_path": None,
        "name": None,
        "acc": None,
        "description": None,
        "date": None,
        "sha512t24u": None,
        "model_length": None,
        "category": None,
        "subcategory": None,
    }
    assert hmms_object.source_counter == 0
    assert hmms_object.subcategory is None
    assert hmms_object.temp_list == []
    assert hmms_object.temp_list2 == []
    assert hmms_object.total_hmms_counter == 0


def test_check():
    hmms_object = HMMParser()
    input_dir = os.path.dirname(pks)
    hmms_object.parse_single_model_file(
        input_dir=input_dir,
        input_path=pks,
        model_source="bigslice",
        input_rel_path="a1/a2",
    )
    assert hmms_object.total_hmms_counter == 26
    assert hmms_object.source_counter == 26
    assert str(hmms_object.rel_path) == "a1/a2"
    hasher = hashlib.sha256()
    hasher.update(str(hmms_object.model_text_dict).encode("utf-8"))
    assert (
        hasher.hexdigest()
        == "bed53994b8be2ffa77d63d3ad2d06ecf9995de53171c52f8ff9e7dacaca0f067"
    )
    hmms_object.change_model_name_to_hash()
    hasher = hashlib.sha256()
    hasher.update(str(hmms_object.model_text_dict).encode("utf-8"))
    assert (
        hasher.hexdigest()
        == "d4b7c8f1a09152f4a248535c118b317c6ffd94c8b9649c8979938e37e2d18745"
    )
    assert (
        hmms_object.model_text_dict[0][1] == "NAME  SMxGr2PCmKBi19aaFwA6UoZHZSUZ5zH9\n"
    )


def test_create_nr_hmm_dict():
    hmms_object = HMMParser()
    input_dir = os.path.dirname(pks)
    hmms_object.parse_single_model_file(
        input_dir=input_dir,
        input_path=pks,
        model_source="bigslice",
        input_rel_path="a1/a2",
    )
    assert hmms_object.nr_models == {}
    hmms_object.create_nr_hmm_dict()
    hasher = hashlib.sha256()
    hasher.update(str(hmms_object.nr_models).encode("utf-8"))
    assert (
        hasher.hexdigest()
        == "97277a445014736eb4d909449003a1b3fcb356f6c450c4719704711527fcc5ad"
    )


def test_bigslice_parse():
    hmms_object = HMMParser()
    input_dir = os.path.dirname(pks)
    hmms_object.parse_single_model_file(
        input_dir=input_dir,
        input_path=pks,
        model_source="bigslice",
        input_rel_path="a1/a2",
    )
    assert hmms_object.model_info_dict[0]["category"] == "a2"
    assert hmms_object.model_info_dict[0]["subcategory"] is None


def test_prism_parse():
    hmms_object = HMMParser()
    input_dir = os.path.dirname(pks)
    hmms_object.parse_single_model_file(
        input_dir=input_dir,
        input_path=pks,
        model_source="prism",
        input_rel_path="a1/a2",
    )
    assert hmms_object.model_info_dict[0]["category"] == "a2"
    assert hmms_object.model_info_dict[0]["subcategory"] is None


def test_antismash_parse():
    def check_antismash(input_rel_path, category, subcategory):
        hmms_object = HMMParser()
        input_dir = os.path.dirname(pks)
        hmms_object.parse_single_model_file(
            input_dir=input_dir,
            input_path=pks,
            model_source="antismash",
            input_rel_path=input_rel_path,
        )
        assert hmms_object.model_info_dict[0]["category"] == category
        if subcategory is None:
            assert hmms_object.model_info_dict[0]["subcategory"] is None
        else:
            assert hmms_object.model_info_dict[0]["subcategory"] == subcategory

    check_antismash(
        input_rel_path="a1/a2/modules/sactipeptides/a1/a2/a3",
        category="sactipeptides",
        subcategory="a2",
    )
    check_antismash(
        input_rel_path="a1/a2/modules/thiopeptides/a1/a2/a3",
        category="thiopeptides",
        subcategory="a2",
    )
    check_antismash(
        input_rel_path="a1/a2/modules/lanthipeptides/a1/a2/a3",
        category="lanthipeptides",
        subcategory="a2",
    )

    check_antismash(
        input_rel_path="a1/a2/detection/hmm_detection/a2/a3",
        category="hmm_detection",
        subcategory=None,
    )
    check_antismash(
        input_rel_path="a1/a2/detection/genefunctions/a2/a3",
        category="genefunctions",
        subcategory=None,
    )
    check_antismash(
        input_rel_path="a1/a2/detection/nrps_pks_domains/a2/dockingdomains.hmm",
        category="nrps_pks_domains",
        subcategory="dockingdomains",
    )
