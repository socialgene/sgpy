import tempfile

from socialgene.base.socialgene import SocialGene


def test_export_all_domains_as_tsv():
    expected = [
        "0hMjYRUCOMiDkJnVKlZ4QVMGhG8mkwdb\t1\t1\t100\t1.1\t0\t-1\t1.1\t1.1\t1.1\t1\t100\t1\t100\tTrue\n",
        "0hMjYRUCOMiDkJnVKlZ4QVMGhG8mkwdb\t10\t1\t100\t1.1\t1.1\t1.0\t1.1\t1.1\t1.1\t1\t100\t1\t100\tFalse\n",
        "0hMjYRUCOMiDkJnVKlZ4QVMGhG8mkwdb\t11\t1\t100\t1.1\t0\t0\t1.1\t1.1\t1.1\t1\t100\t1\t100\tTrue\n",
        "0hMjYRUCOMiDkJnVKlZ4QVMGhG8mkwdb\t12\t1\t100\t1.1\t1.1\t0.0\t1.1\t1.1\t1.1\t1\t100\t1\t100\tFalse\n",
        "0hMjYRUCOMiDkJnVKlZ4QVMGhG8mkwdb\t2\t1\t100\t1.1\t1.1\t0.1\t1.1\t1.1\t1.1\t1\t100\t1\t100\tFalse\n",
        "0hMjYRUCOMiDkJnVKlZ4QVMGhG8mkwdb\t3\t1\t100\t1.1\t0\t2\t1.1\t1.1\t1.1\t1\t100\t1\t100\tTrue\n",
        "0hMjYRUCOMiDkJnVKlZ4QVMGhG8mkwdb\t4\t1\t100\t1.1\t1.1\t100.0\t1.1\t1.1\t1.1\t1\t100\t1\t100\tFalse\n",
        "0hMjYRUCOMiDkJnVKlZ4QVMGhG8mkwdb\t5\t1\t100\t1.1\t0\t3\t1.1\t1.1\t1.1\t1\t100\t1\t100\tTrue\n",
        "0hMjYRUCOMiDkJnVKlZ4QVMGhG8mkwdb\t6\t1\t100\t1.1\t1.1\t1000.0\t1.1\t1.1\t1.1\t1\t100\t1\t100\tFalse\n",
        "0hMjYRUCOMiDkJnVKlZ4QVMGhG8mkwdb\t7\t1\t100\t1.1\t0\t0\t1.1\t1.1\t1.1\t1\t100\t1\t100\tTrue\n",
        "0hMjYRUCOMiDkJnVKlZ4QVMGhG8mkwdb\t8\t1\t100\t1.1\t1.1\t2.0\t1.1\t1.1\t1.1\t1\t100\t1\t100\tFalse\n",
        "0hMjYRUCOMiDkJnVKlZ4QVMGhG8mkwdb\t9\t1\t100\t1.1\t0\t0\t1.1\t1.1\t1.1\t1\t100\t1\t100\tTrue\n",
    ]
    sg_obj = SocialGene()
    sg_obj.add_protein(
        sequence="ARNDCQEGHILKMFPSTWYVXZJU",
        description="description",
        external_protein_id="external_protein_id",
    )
    base_dict = {
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
        {"hmm_id": "1", "exponentialized": True, "i_evalue": 0.1},
        {"hmm_id": "2", "exponentialized": False, "i_evalue": 0.1},
        {"hmm_id": "3", "exponentialized": True, "i_evalue": 100},
        {"hmm_id": "4", "exponentialized": False, "i_evalue": 100},
        {
            "hmm_id": "5",
            "exponentialized": True,
            "i_evalue": 1000,
        },  # should fail to pass threshold
        {
            "hmm_id": "6",
            "exponentialized": False,
            "i_evalue": 1000,
        },  # should fail to pass threshold
        {"hmm_id": "7", "exponentialized": True, "i_evalue": 2},
        {"hmm_id": "8", "exponentialized": False, "i_evalue": 2},
        {"hmm_id": "9", "exponentialized": True, "i_evalue": 1},
        {"hmm_id": "10", "exponentialized": False, "i_evalue": 1},
        {"hmm_id": "11", "exponentialized": True, "i_evalue": 0},
        {"hmm_id": "12", "exponentialized": False, "i_evalue": 0},
    )
    for i in [base_dict | i for i in tempd]:
        sg_obj.proteins["0hMjYRUCOMiDkJnVKlZ4QVMGhG8mkwdb"].add_domain(**i)

    with tempfile.NamedTemporaryFile() as temp_path:
        sg_obj.export_all_domains_as_tsv(temp_path.name)
        with open(temp_path.name) as h:
            z = h.readlines()
    # domains are stored internally as sets, so may have different order in tsv
    assert sorted(z) == sorted(expected)


TODO: 1
