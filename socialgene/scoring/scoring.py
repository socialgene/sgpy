from collections import OrderedDict

from textdistance import jaccard, levenshtein


def mod_score(input_list_1, input_list_2):
    """Combine Edit and Composition Distances (Levenshtein and Jaccard)
    Args:
        input_list_1 (List): list to compare to input_list_2
        input_list_2 (List): list to compare to input_list_1
    Returns:
        dict: {l1, l2, levenshtein, jaccard, mod_score}; mod_score -> 2 = Perfectly similar; otherwise (1/Levenshtein + Jaccard)
    """
    # If either protein contains no HMM annotations,
    # return a mod score with the worst scores possible
    if not isinstance(input_list_1, list) or not isinstance(input_list_2, list):
        raise TypeError()
    # have to take length here because it might be artificially changed later
    length_input_list_1 = len(input_list_1)
    length_input_list_2 = len(input_list_2)
    if not input_list_1 or not input_list_2:
        # If either protein contains no HMM annotations,
        # return a mod score with the worst scores possible
        input_list_1 = [0]
        input_list_2 = [1]
    # levenshtein_score: 0=most similar; increases based on number of edits
    levenshtein_score = levenshtein(input_list_1, input_list_2)
    # convert int64 to int
    levenshtein_score = int(levenshtein_score)
    # jaccard_score: 1= most similar; 0=most dissimilar
    jaccard_score = jaccard(input_list_1, input_list_2)
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

    return OrderedDict(
        {
            "l1": length_input_list_1,
            "l2": length_input_list_2,
            "levenshtein": round(mod_levenshtein_score, 2),
            "jaccard": round(jaccard_score, 2),
            "mod_score": round(mod_score_value, 2),
        }
    )
