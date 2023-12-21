from collections import namedtuple

from textdistance import jaccard, levenshtein

from socialgene.base.molbio import Protein

# Note: there is technically no query/target as the comparison methods are order-independent
# this is just to keep things consistent between tools
_mod_score_tupler = namedtuple(
    "protein_comparison_modscore",
    (
        "query",
        "target",
        "query_n_domains",
        "target_n_domains",
        "levenshtein",
        "jaccard",
        "mod_score",
    ),
)


def mod_score(p1, p2):
    """Combine Edit and Composition Distances (Levenshtein and Jaccard)
    Important: If protein sequence hashes are identical then a perfect score is returned without comparing domains.
    Args:
        p1 (Protein): SocialGene Protein Class
        p1 (Protein): SocialGene Protein Class
    Returns:
        dict: {ProteinClass, ProteinClass, levenshtein, jaccard, mod_score}; mod_score -> 2 = Perfectly similar; otherwise (1/Levenshtein + Jaccard)
    """
    # If either protein contains no HMM annotations,
    # return a mod score with the worst scores possible
    if not isinstance(p1, Protein) and not isinstance(p2, Protein):
        raise TypeError(f"p1 type: {type(p1)}; p2 type: {type(p2)}")
    length_input_list_1 = len(p1.domains)
    length_input_list_2 = len(p2.domains)
    if p1 == p2:
        return _mod_score_tupler(
            p1,
            p2,
            length_input_list_1,
            length_input_list_2,
            round(0, 2),
            round(1, 2),
            round(1.5, 2),
        )
    if not p1.domains or not p2.domains:
        # If either protein contains no HMM annotations,
        # return a mod score with the worst scores possible
        return _mod_score_tupler(
            p1,
            p2,
            length_input_list_1,
            length_input_list_2,
            round(100, 2),
            round(0, 2),
            round(0, 2),
        )
    # levenshtein_score: 0=most similar; increases based on number of edits
    levenshtein_score = levenshtein(p1.domain_vector, p2.domain_vector)
    # convert int64 to int
    levenshtein_score = int(levenshtein_score)
    # jaccard_score: 1= most similar; 0=most dissimilar
    jaccard_score = jaccard(set(p1.domain_vector), set(p2.domain_vector))
    # change levenshtein so 0=worst, 1 = best
    max_edit = max(length_input_list_1, length_input_list_2)
    if max_edit == 0:
        # Neither protein had HMMs
        mod_score_value = 0
        mod_levenshtein_score = 0
    elif length_input_list_1 == 0 or length_input_list_2 == 0:
        # One protein didn't have HMMs
        mod_score_value = 0
        mod_levenshtein_score = 0
    else:
        mod_levenshtein_score = 1 - (levenshtein_score / max_edit)
        if jaccard_score == 0:
            mod_score_value = 0
        else:
            mod_score_value = (jaccard_score * 0.5) + mod_levenshtein_score

    return _mod_score_tupler(
        p1,
        p2,
        length_input_list_1,
        length_input_list_2,
        round(mod_levenshtein_score, 2),
        round(jaccard_score, 2),
        round(mod_score_value, 2),
    )
