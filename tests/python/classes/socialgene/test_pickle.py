import tempfile
from socialgene.base.socialgene import SocialGene
from socialgene.config import env_vars


def test_read_and_writing_pickled_sgobject():
    sg_obj = SocialGene()
    env_vars["HMMSEARCH_IEVALUE"] = 1000
    _ = sg_obj.add_protein(
        sequence="ARNDCQEGHILKMFPSTWYVXZJU",
        description="description",
        other_id="other_id",
    )
    base_dict = {
        "hmm_id": "hmm_id",
        "env_from": 1,
        "env_to": 100,
        "seq_pro_score": 1.1,
        "evalue": 1.1,
        "domain_bias": 1.1,
        "domain_score": 1.1,
        "seq_pro_bias": 1.1,
        "hmm_from": 1,
        "hmm_to": 100,
        "ali_from": 1,
        "ali_to": 100,
    }
    tempd = (
        {"exponentialized": True, "i_evalue": 0.1},
        {"exponentialized": False, "i_evalue": 2},
        {"exponentialized": False, "i_evalue": 1},
        {"exponentialized": False, "i_evalue": 0},
        {"exponentialized": False, "i_evalue": 0.1},
    )
    for i in [base_dict | i for i in tempd]:
        sg_obj.proteins["0hMjYRUCOMiDkJnVKlZ4QVMGhG8mkwdb"].add_domain(**i)

    sg_obj.add_assembly("myassembly")
    sg_obj.assemblies["myassembly"].add_locus(
        "my_locus",
    )
    sg_obj.assemblies["myassembly"].loci["my_locus"].add_feature(
        id="feature_id1",
        type="protein",
        start=1,
        end=10,
        strand=1,
    )
    sg_obj.assemblies["myassembly"].loci["my_locus"].add_feature(
        id="feature_id2",
        type="not_a_prot",
        start=1,
        end=10,
        strand=1,
    )
    temp_path = tempfile.NamedTemporaryFile()
    # write the pickled sg object
    sg_obj.ferment_pickle(temp_path.name)
    # Test reading the pickle back in
    new_sg_obj = SocialGene.eat_pickle(temp_path.name)
    assert (
        new_sg_obj.proteins["0hMjYRUCOMiDkJnVKlZ4QVMGhG8mkwdb"].other_id == "other_id"
    )
    assert (
        new_sg_obj.proteins["0hMjYRUCOMiDkJnVKlZ4QVMGhG8mkwdb"].sequence
        == "ARNDCQEGHILKMFPSTWYVXZJU"
    )
    assert (
        new_sg_obj.proteins["0hMjYRUCOMiDkJnVKlZ4QVMGhG8mkwdb"].description
        == "description"
    )
    assert (
        new_sg_obj.proteins["0hMjYRUCOMiDkJnVKlZ4QVMGhG8mkwdb"].hash_id
        == "0hMjYRUCOMiDkJnVKlZ4QVMGhG8mkwdb"
    )
    assert sorted(
        [i.id for i in new_sg_obj.assemblies["myassembly"].loci["my_locus"].features]
    ) == ["feature_id1", "feature_id2"]
