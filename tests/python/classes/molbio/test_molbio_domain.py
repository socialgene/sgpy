import pytest
from socialgene.base.molbio import Domain
from socialgene.config import env_vars


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


temp = (
    [{"exponentialized": False, "i_evalue": 1000}, False],
    [{"exponentialized": False, "i_evalue": 10}, False],
    [{"exponentialized": False, "i_evalue": 1.1}, False],
    [{"exponentialized": False, "i_evalue": 1}, False],
    [{"exponentialized": False, "i_evalue": 0.1}, False],
    [{"exponentialized": False, "i_evalue": 0}, False],
)

combos = [(base_dict | i[0], i[1]) for i in temp]
del temp
del base_dict


@pytest.mark.parametrize("a,expected", combos)
def test_domain_1(a, expected):
    env_vars["HMMSEARCH_IEVALUE"] = 0.2
    temp = Domain(**a)
    assert temp.domain_within_threshold() == expected


#     assert temp.get_dict() == OrderedDict(
#         [
#             ("hmm_id", "hmm_id"),
#             ("env_from", 1),
#             ("env_to", 100),
#             ("seq_pro_score", 1.1),
#             ("evalue", 0),
#             ("i_evalue", 0),
#             ("domain_bias", 1.1),
#             ("domain_score", 1.1),
#             ("seq_pro_bias", 1.1),
#             ("hmm_from", 1),
#             ("hmm_to", 100),
#             ("ali_from", 1),
#             ("ali_to", 100),
#         ]
#     )
#     env_vars["HMMSEARCH_IEVALUE"] = 0
#     assert temp.domain_within_threshold() is True
#     env_vars["HMMSEARCH_IEVALUE"] = 0.1
#     assert temp.domain_within_threshold() is False
#     env_vars["HMMSEARCH_IEVALUE"] = 1
#     assert temp.domain_within_threshold() is True
#     env_vars["HMMSEARCH_IEVALUE"] = 0.0001
#     assert temp.domain_within_threshold() is False


# def test_domain_2_fail():
#     with pytest.raises(ValueError):
#         env_vars["HMMSEARCH_IEVALUE"] = -1
#         temp.domain_within_threshold()


# def test_exponentialized():
#     temp = Domain(
#         hmm_id="hmm_id",
#         exponentialized=False,
#         env_from=1,
#         env_to=100,
#         seq_pro_score=1.1,
#         evalue=1.1,
#         i_evalue=1.1,
#         domain_bias=1.1,
#         domain_score=1.1,
#         seq_pro_bias=1.1,
#         hmm_from=1,
#         hmm_to=100,
#         ali_from=1,
#         ali_to=100,
#     )
#     assert temp.get_dict() == OrderedDict(
#         [
#             ("hmm_id", "hmm_id"),
#             ("env_from", 1),
#             ("env_to", 100),
#             ("seq_pro_score", 1.1),
#             ("evalue", 0),
#             ("i_evalue", 0),
#             ("domain_bias", 1.1),
#             ("domain_score", 1.1),
#             ("seq_pro_bias", 1.1),
#             ("hmm_from", 1),
#             ("hmm_to", 100),
#             ("ali_from", 1),
#             ("ali_to", 100),
#         ]
#     )
#     temp = Domain(
#         hmm_id="hmm_id",
#         exponentialized=True,
#         env_from=1,
#         env_to=100,
#         seq_pro_score=1.1,
#         evalue=1.1,
#         i_evalue=1.1,
#         domain_bias=1.1,
#         domain_score=1.1,
#         seq_pro_bias=1.1,
#         hmm_from=1,
#         hmm_to=100,
#         ali_from=1,
#         ali_to=100,
#     )
#     assert temp.get_dict() == OrderedDict(
#         [
#             ("hmm_id", "hmm_id"),
#             ("env_from", 1),
#             ("env_to", 100),
#             ("seq_pro_score", 1.1),
#             ("evalue", 1),
#             ("i_evalue", 1),
#             ("domain_bias", 1.1),
#             ("domain_score", 1.1),
#             ("seq_pro_bias", 1.1),
#             ("hmm_from", 1),
#             ("hmm_to", 100),
#             ("ali_from", 1),
#             ("ali_to", 100),
#         ]
#     )
