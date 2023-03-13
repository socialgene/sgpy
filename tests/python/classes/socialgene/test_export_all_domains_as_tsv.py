import tempfile
from socialgene.base.socialgene import SocialGene
from socialgene.config import env_vars


def test_export_all_domains_as_tsv():
    expected = [
        "0hMjYRUCOMiDkJnVKlZ4QVMGhG8mkwdb\thmm_id\t1\t100\t1.1\t1\t2\t1.1\t1.1\t1.1\t1\t100\t1\t100\n",
        "0hMjYRUCOMiDkJnVKlZ4QVMGhG8mkwdb\thmm_id\t1\t100\t1.1\t0\t-1\t1.1\t1.1\t1.1\t1\t100\t1\t100\n",
        "0hMjYRUCOMiDkJnVKlZ4QVMGhG8mkwdb\thmm_id\t1\t100\t1.1\t1\t0\t1.1\t1.1\t1.1\t1\t100\t1\t100\n",
        "0hMjYRUCOMiDkJnVKlZ4QVMGhG8mkwdb\thmm_id\t1\t100\t1.1\t1\t1\t1.1\t1.1\t1.1\t1\t100\t1\t100\n",
    ]
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

    with tempfile.NamedTemporaryFile() as temp_path:
        sg_obj.export_all_domains_as_tsv(temp_path.name)
        with open(temp_path.name) as h:
            z = h.readlines()
    # domains are stored internally as sets, so may have different order in tsv
    assert sorted(z) == sorted(expected)
