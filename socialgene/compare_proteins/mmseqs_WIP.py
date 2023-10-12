# flake8: noqa
from itertools import combinations, product
from multiprocessing import Pool
from pathlib import Path

import pandas as pd

from socialgene.compare_proteins.base_class import CompareProteinsBaseClass
from socialgene.mmseqs.create_database import create_database
from socialgene.mmseqs.search import search
from socialgene.mmseqs.subset_database import createsubdb

target_proteins = list(sg.proteins.values())

# with tempfile.TemporaryDirectory() as tmpdirname:

tmpdirname = "/home/chase/Downloads/work"
fasta_path = Path(tmpdirname, "target_fasta.fa")
mmseqs_path = "/home/chase/Downloads/work/mmseqsdb"

with open(fasta_path, "w") as handle:
    handle.writelines(
        (i.fasta_string_defline_hash_id for i in target_proteins if i.sequence)
    )


create_database(fasta_path, mmseqs_path)

mmseqs_lookup = "/home/chase/Downloads/work/mmseqsdb.lookup"
mmseqs_output_subset_ids = "/home/chase/Downloads/work/subset_ids"
mmseqs_subset_db = "/home/chase/Downloads/work/subset_db"

# create a subdb using protein hash ids
protids = [i.hash_id for i in target_proteins][0:4]
with open(mmseqs_lookup, "r") as mml:
    with open(mmseqs_output_subset_ids, "w") as h:
        for i in mml:
            temp = i.split("\t")
            if temp[1] in protids:
                z = "\t"
                h.writelines(f"{i.split(z)[0]}\n")

createsubdb(olddb=mmseqs_path, newdb=mmseqs_subset_db, idfile=mmseqs_output_subset_ids)


a = search(fasta_path, mmseqs_subset_db)


class CompareDiamond(CompareProteinsBaseClass):
    def __init__(self):
        self.protein_comparisons = set()

    @property
    def mod_score_df(self):
        return (
            pd.DataFrame(self.protein_comparisons)
            .sort_values(by=["mod_score"], ascending=False)
            .reset_index(inplace=False, drop=True)
        )

    @property
    def df(self):
        return self.mod_score_df.filter(["query", "target", "mod_score"]).rename(
            columns={"mod_score": "score"}
        )

    def compare_one_to_one(self, p1, p2):
        return self.calculate_mod_score(p1, p2)

    def compare_one_to_many(self, p1_obj, p2_obj_list):
        for i in p2_obj_list:
            temp = self.calculate_mod_score(p1_obj, i)
            if temp.jaccard > 0:
                self.protein_comparisons.add(temp)

    def compare_many_to_many(self, p1_obj_list, p2_obj_list):
        for i1, i2 in product(p1_obj_list, p2_obj_list):
            temp = self.calculate_mod_score(i1, i2)
            if temp.jaccard > 0:
                self.protein_comparisons.add(temp)

    def compare_all_to_all(self, p1_obj_list):
        for i1, i2 in combinations(p1_obj_list, 2):
            temp = self.calculate_mod_score(i1, i2)
            if temp.jaccard > 0:
                temp
                self.protein_comparisons.add(temp)

    def compare_all_to_all_parallel(self, p1_obj_list, cpus=1):
        # have to use _calculate_mod_score_not_named because named tuple can't pickle "protein_comparison_modscore"
        with Pool(cpus) as p:
            for i in p.starmap(
                _calculate_mod_score_not_named,
                combinations(p1_obj_list, 2),
            ):
                # i[5] == jaccard
                if i[5] > 0.001:
                    self.protein_comparisons.add(_create_tuple(*i))
