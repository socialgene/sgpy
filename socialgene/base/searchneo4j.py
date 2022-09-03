#!/usr/bin/env python

# python dependencies
from collections import defaultdict

# external dependencies

# internal dependencies
from socialgene.base.socialgene import SocialGene
from socialgene.base.compare_protein import CompareProtein
from socialgene.neo4j.single_protein_search import SingleProteinSearch


class SearchNeo4j(CompareProtein, SingleProteinSearch):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # used for parsing/holding input protein(s)
        self.input_sg_object = SocialGene()
        # used for holding result protein(s)
        self.result_sg_object = SocialGene()
        # used for holding query/match dict {input_protein1: [result_protein1, result_protein2]}
        self.query_and_match = defaultdict(set)
        self.query_subject_levenshtein_jaccard_mod_score_table = None
