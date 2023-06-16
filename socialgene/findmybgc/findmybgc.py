from collections import defaultdict

from socialgene.base.molbio import Protein
from socialgene.base.socialgene import SocialGene
from socialgene.neo4j.neo4j import Neo4jQuery
from socialgene.utils.logging import log


class SingleProteinSearch(SocialGene):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.query_targets = dict()

    def search_a_single_protein(
        self, protein_object: Protein, modscore_threshold: float = 0.5, max_matches=None
    ):
        """Search a database for similar proteins (crude search using domain content only)

        Args:
            protein_object (Protein): socialgene protein object
            modscore_threshold (float, optional): will ignore comparisons with a mod_score below this. Defaults to 0.5.
            max_matches (_type_, optional): Limit return to x-matches. Defaults to None.
        """

        log.info(f"Input {protein_object.external_protein_id}: begin database search")
        self.query_targets[protein_object.hash_id] = set()
        targets = [
            i["target"]
            for i in Neo4jQuery.query_neo4j(
                cypher_name="search_a_single_protein",
                param=protein_object.get_domain_vector(only_unique=True),
            )
        ]
        if not targets:
            log.warning(f"Input {protein_object.external_protein_id}: no matches found")
            return
        log.info(
            f"Input {protein_object.external_protein_id}: first pass database search returned {len(targets)} proteins."
        )
        if max_matches and len(targets) > max_matches:
            # TODO: this should be done before this by limiting the neo4j results
            log.warning(
                f"Input {protein_object.external_protein_id}: Truncating result list to {max_matches} matches"
            )
            targets = targets[0:max_matches]
        for hash_id in targets:
            _ = self.add_protein(hash_id=hash_id)
        self.get_protein_domains_from_db(protein_id_list=targets)
        for target in targets:
            score = self._calculate_mod_score_from_protein_class(
                protein_object, self.proteins[target]
            )
            if score[6] > modscore_threshold:
                self.protein_comparison.append(score)
                self.query_targets[protein_object.hash_id].add(target)
            else:
                del self.proteins[target]
        for i in self.proteins.values():
            i.only_id()
        self.fill_locus_assembly_from_db()
        log.info(
            f"Input {protein_object.external_protein_id}: second pass returned {len(self.proteins)} proteins and {len(self.assemblies)} assemblies"
        )


class NewFindMyBGC:
    def __init__(self):
        self.input_sg_object = SocialGene()
        self.big_sg_list = dict()
        self.count_of_inputs_per_assembly = defaultdict(int)
        self.count_of_inputs_per_locus = {}
        self.query_targets = dict()

    def parse_input(self, input_path):
        """Use SocialGene() class's parser on the input BGC
        Args:
            input_path (str): input BGC filepath
        """
        self.input_sg_object.parse(input_path)

    def annotate(self, **kwargs):
        """Annotate the input proteins with HMMs"""
        self.input_sg_object.annotate(**kwargs)

    def search_database(self, **kwargs):
        """Search the database, one protein at a time"""
        for protein_k, protein_v in self.input_sg_object.proteins.items():
            self.big_sg_list[protein_k] = SingleProteinSearch()
            self.big_sg_list[protein_k].search_a_single_protein(protein_v, **kwargs)

    def count_inputs_per_assembly_before_merge(self):
        self.count_of_inputs_per_assembly = defaultdict(int)
        for key in self.big_sg_list.values():
            for i in key.assemblies.keys():
                self.count_of_inputs_per_assembly[i] += 1

    def count_inputs_per_assembly_after_merge(self):
        self.count_of_inputs_per_assembly = defaultdict(int)
        a = dict()
        queries = list(self.query_targets.keys())
        for i in queries:
            a[i] = set()
        for k1, v1 in self.query_targets.items():
            for k2, v2 in self.input_sg_object.assemblies.items():
                for k3, v3 in v2.loci.items():
                    if [i.protein_hash for i in v3.features if i.protein_hash in v1]:
                        a[k1].add(k2)
        for k, v in a.items():
            for i in v:
                self.count_of_inputs_per_assembly[i] += 1
        return self.count_of_inputs_per_assembly

    def count_inputs_per_locus(self):
        self.count_of_inputs_per_locus = {}
        for v1 in self.big_sg_list.values():
            for k2, v2 in v1.assemblies.items():
                if k2 not in self.count_of_inputs_per_locus:
                    self.count_of_inputs_per_locus[k2] = {}
                for k3 in v2.loci.keys():
                    if k3 not in self.count_of_inputs_per_locus[k2]:
                        self.count_of_inputs_per_locus[k2][k3] = 0
                    self.count_of_inputs_per_locus[k2][k3] += 1

    def drop_assemblies(self, min_matches):
        """Drop assemblies based on the number of input proteins matched to an assembly

        Args:
            min_matches (int): Minimum number of input proteins that must have a match in an assembly
        """
        # used for dropping assemblies with few matches to the input BGC
        for i in [
            k for k, v in self.count_of_inputs_per_assembly.items() if v < min_matches
        ]:
            for v in self.big_sg_list.values():
                if i in v.assemblies:
                    del v.assemblies[i]
        for i in [k for k, v in self.big_sg_list.items() if not v]:
            # delete any inputs in big_sg_list if there are no assemblies left
            del self.big_sg_list[i]

    def drop_assemblies_by_count22222222(self, min_matches):
        """Drop assemblies based on the number of input proteins matched to an assembly

        Args:
            min_matches (int): Minimum number of input proteins that must have a match in an assembly
        """
        # used for dropping assemblies with few matches to the input BGC
        for i in [
            k for k, v in self.count_of_inputs_per_assembly.items() if v < min_matches
        ]:
            if i in self.input_sg_object.assemblies:
                del self.input_sg_object.assemblies[i]

    def _merge_all(self):
        for k, v in self.big_sg_list.items():
            if isinstance(v, SingleProteinSearch):
                self.input_sg_object + v
                self.query_targets.update(self.big_sg_list[k].query_targets)
            # self.big_sg_list[k] = None

    def calculate_mod_scores(self):
        for query, targets in self.query_targets.items():
            for target in targets:
                self.input_sg_object.protein_comparison.append(
                    self.input_sg_object._calculate_mod_score_from_protein_class(
                        protein_1=self.input_sg_object.proteins[query],
                        protein_2=self.input_sg_object.proteins[target],
                    )
                )


# zz = NewFindMyBGC()
# zz.parse(gbk_path)
# zz.annotate(
#     hmm_filepath=hmm_filepath,
#     cpus=1,
#     use_neo4j_precalc=True,
# )
# zz.loopit()
# zz.big_sg_list

# zz.count_inputs_per_assembly()
# zz.big_sg_list
# zz.count_inputs_per_locus()
# zz.big_sg_list
# zz.count_of_inputs_per_assembly
# # zz.drop_assemblies_by_count(17)
# zz.count_of_inputs_per_locus
# zz.intactness()
# zz.big_sg_list

# zz.input_sg_object.assemblies

# zz._merge_all()
# zz.break_loci()


# def only_largest_bgc(self):
#     log.info(
#         "Retaining single locus per assembly that contains the largest set of genes, this is a destructive action"
#     )
#     starting_loci_count = sum(
#         [len(i.loci) for i in self.result_sg_object.assemblies.values()]
#     )
#     # TODO: don't copy
#     result_holder = deepcopy(self.result_sg_object.assemblies)
#     # loop through assemblies, only keep loci with most features
#     for k, v in self.result_sg_object.assemblies.items():
#         # get number of features of all loc in a single assembly
#         a_dictionary = {k: len(v.features) for k, v in v.loci.items()}
#         # find locus with most features
#         max_key = max(a_dictionary, key=a_dictionary.get)
#         # remove all other loci that don't have the max features
#         for k2, v2 in v.loci.items():
#             if k2 != max_key:
#                 del result_holder[k].loci[k2]
#     self.result_sg_object.assemblies = result_holder
#     del result_holder
#     ending_loci_count = sum(
#         [len(i.loci) for i in self.result_sg_object.assemblies.values()]
#     )
#     log.info(f"Removed {starting_loci_count - ending_loci_count} loci ")

# def min_genes(self, min_genes):
#     n_start = len(self.result_sg_object.assemblies)
#     # TODO: don't copy
#     result_holder = deepcopy(self.result_sg_object.assemblies)
#     # remove assemblies if they don't have a locus with >= min_genes
#     for k, v in self.result_sg_object.assemblies.items():
#         n_genes = 0
#         for v2 in v.loci.values():
#             loc_len = len(v2.features)
#             if loc_len > n_genes:
#                 n_genes = loc_len
#         if n_genes < min_genes:
#             temp = result_holder.pop(k)
#     self.result_sg_object.assemblies = result_holder
#     del result_holder
#     n_end = len(self.result_sg_object.assemblies)
#     log.info(f"Removed:{n_start-n_end} assemblies; {n_end} assemblies remain")

# def calculate_mod_scores(self):

#     for query, targets in self.found_proteins.items():
#         for target in targets:
#             self.protein_comparison.append(
#                 self._calculate_mod_score_from_protein_class(
#                     protein_1=self.input_sg_object.proteins[query],
#                     protein_2=self.result_sg_object.proteins[target],
#                 )
#             )
