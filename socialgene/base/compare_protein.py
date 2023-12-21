from itertools import combinations_with_replacement, product
from multiprocessing import Pool, cpu_count

import pandas as pd

from socialgene.compare_proteins.hmm_scoring import mod_score
from socialgene.neo4j.neo4j import Neo4jQuery
from socialgene.utils.logging import log


def _mod_return(i1, i2):
    """for running mod_score() in parallel"""
    return (
        i1[0],  # hash of protein 1
        i2[0],  # hash of protein 2
        *(
            mod_score(
                i1[1],
                i2[1].domain_vector,
            ).values()
        ),
    )


def append_or_not(result_list, result, append=False):
    if append:
        result_list.extend(result)
    else:
        result_list = result


class CompareProtein(Neo4jQuery):
    def __init__(
        self,
    ):
        super().__init__()
        self.protein_comparison = []

    def _calculate_mod_score_from_protein_class(
        self,
        protein_1=None,
        protein_2=None,
        **kwargs,
    ):
        """Compare two proteins, using protein class as input

        Args:
            protein_1 (Protein, optional): socialgene Protein object. Defaults to None.
            protein_2 (Protein, optional): socialgene Protein object. Defaults to None.
            append (bool, optional): Should results be appended to previously calculated comparisons? Defaults to False.
            verbose (bool, optional): If `True` then additional info will be printed. Defaults to False.
        """

        return self._calculate_mod_score_from_domain_lists(
            protein_id_1=protein_1.uid,
            protein_id_2=protein_2.uid,
            input_list_1=protein_1.domain_vector,
            input_list_2=protein_2.domain_vector,
            **kwargs,
        )

    @staticmethod
    def _calculate_mod_score_from_domain_lists(
        protein_id_1,
        protein_id_2,
        input_list_1,
        input_list_2,
    ):
        return [
            protein_id_1,
            protein_id_2,
            *(mod_score(input_list_1=input_list_1, input_list_2=input_list_2).values()),
        ]

    @staticmethod
    def compare_two_proteins(protein_1, protein_2):
        """
        dict: {l1, l2, levenshtein, jaccard, mod_score}; mod_score -> 2 = Perfectly similar; otherwise (1/Levenshtein + Jaccard)
        Returns:

          protein_2: An instance of a protein object representing the second protein.
          protein_1: An instance of a protein object representing the first protein.
        Args:

        `mod_score` function.
        The function `compare_two_proteins` compares the domain vectors of two proteins using the
        """

        return mod_score(protein_1.domain_vector, protein_2.domain_vector)

    def bro(self, queries, targets, **kwargs):
        """Compare a list fo proteins to another list of proteins

        Args:
            query (List(str)): list of protein hashes
            targets (List(str): list of protein hashes
        """
        self.compare_proteins(product(queries, targets), **kwargs)

    def compare_proteins(
        self,
        uid_list_of_tuples=None,
        cpus=None,
        append=False,
        verbose=False,
    ):
        """Calculate similarity between proteins based on their domain content

        Args:
            uid_list_of_tuples (str, optional): list of two-ples of protein hashes; use `None` to calculate all-vs-all. Defaults to None.
            cpus (int, optional): If >1 then multiprocessing may/may-not be used to speed things up. Defaults to None.
            append (bool, optional): Should results be appended to previously calculated comparisons? Defaults to False.
            verbose (bool, optional): If `True` then additional info will be printed. Defaults to False.
        """
        if uid_list_of_tuples is not None:
            for id_pair in uid_list_of_tuples:
                _temp = self._calculate_mod_score_from_domain_lists(
                    protein_id_1=id_pair[0],
                    protein_id_2=id_pair[1],
                    input_list_1=self.proteins[id_pair[0]].domain_vector,
                    input_list_2=self.proteins[id_pair[1]].domain_vector,
                )
                if append:
                    self.protein_comparison.append(_temp)
                else:
                    return _temp

        elif uid_list_of_tuples is None:
            if len(self.proteins) < 10:
                if verbose:
                    log.info(
                        f"Running all-vs-all comparison of {len(self.proteins)} proteins (list comprehension method)"
                    )
                append_or_not(
                    result_list=self.protein_comparison,
                    result=[
                        _mod_return(i1=i1, i2=i2)
                        for i1, i2 in combinations_with_replacement(
                            self.proteins.items(), r=2
                        )
                    ],
                    append=append,
                )
            else:
                if cpus is None:
                    cpus = cpu_count() - 1
                    if cpus < 1:
                        cpus = 1
                if verbose:
                    log.info(
                        f"Running all-vs-all comparison of {len(self.proteins)} proteins (multiprocessing method, using {cpus}-cores)"
                    )
                with Pool(cpus) as p:
                    append_or_not(
                        result_list=self.protein_comparison,
                        result=p.starmap(
                            _mod_return,
                            combinations_with_replacement(self.proteins.items(), r=2),
                        ),
                        append=append,
                    )

    def protein_comparison_to_df(self):
        if isinstance(self.protein_comparison, list):
            log.info(
                "Converting protein_comparison list to dataframe and sorting by mod_score"
            )
            self.protein_comparison = pd.DataFrame(
                self.protein_comparison,
                columns=[
                    "query",
                    "target",
                    "query_n_domains",
                    "target_n_domains",
                    "levenshtein",
                    "jaccard",
                    "mod_score",
                ],
            )
            self.protein_comparison.sort_values(
                by=["mod_score"], inplace=True, ascending=False
            )
        elif isinstance(self.protein_comparison, pd.DataFrame):
            log.info(
                "protein_comparison is already a pandas dataframe, skipping conversion"
            )
        else:
            log.info("confused.....protein_comparison isn't list or dataframe")

    def protein_comparison_append_species_and_assemblies(self):
        """Take the pandas DF target column, and search neo4j for those protein ids, return each
        assemblies and species name for each protein. Merge it into the pandas DF.
        """
        result = self.query_neo4j(
            cypher_name="get_species_and_assembly_from_protein_list",
            param=list(set(self.protein_comparison["target"])),
        )
        result = pd.DataFrame(result)

        self.protein_comparison = pd.merge(
            self.protein_comparison,
            result,
            how="left",
            on=None,
            left_on="target",
            right_on="protein",
            left_index=False,
            right_index=False,
            sort=False,
            suffixes=("_x", "_y"),
            copy=True,
            indicator=False,
            validate=None,
        )
        self.protein_comparison.drop(
            "protein",
            axis=1,
            inplace=True,
        )

    def protein_comparison_write_json(self, outpath=None, orient="records"):
        # orient=values is smaller output than orient=records
        if outpath is None:
            return self.protein_comparison.to_json(
                path_or_buf=outpath, orient=orient, double_precision=2, force_ascii=True
            )
        else:
            self.protein_comparison.to_json(
                path_or_buf=outpath, orient=orient, double_precision=2, force_ascii=True
            )

    def protein_comparison_sort(self):
        self.protein_comparison.sort_values(
            by=["mod_score"], inplace=True, ascending=False
        )

    def protein_comparison_append_protein_names(self):
        """Take the pandas DF target column, and search neo4j for those protein ids, return each
        assemblies and species name for each protein. Merge it into the pandas DF.
        """
        result = self.query_neo4j(
            cypher_name="get_protein_name",
            param=list(set(self.protein_comparison["target"])),
        )
        result = pd.DataFrame(result)
        self.protein_comparison = pd.merge(
            self.protein_comparison,
            result,
            how="left",
            on=None,
            left_on="target",
            right_on="external_id",
            left_index=False,
            right_index=False,
            sort=False,
            suffixes=("_x", "_y"),
            copy=True,
            indicator=False,
            validate=None,
        )
        self.protein_comparison.drop(
            "external_id",
            axis=1,
            inplace=True,
        )
