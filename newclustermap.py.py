from socialgene.clustermap.clustermap import Clustermap
from socialgene.hmm.hmmer import HMMER
from socialgene.base.socialgene import SocialGene
from pathlib import Path
import pandas as pd
from rich import inspect


sg_object = SocialGene()

sg_object.parse("/home/chase/Documents/data/mibig/3_1/mibig_gbk_3.1/BGC0000001.gbk")
sg_object.parse("/home/chase/Documents/data/mibig/3_1/mibig_gbk_3.1/BGC0000002.gbk")

self = sg_object


class TabularRow:
    # http://useast.ensembl.org/info/website/upload/gff3.html
    __slots__ = [
        "assembly",
        "seqid",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attribute",
    ]

    def __init__(
        self,
        assembly=None,
        seqid=None,
        source=None,
        type=None,
        start=None,
        end=None,
        score=None,
        strand=None,
        phase=None,
        attribute=None,
    ) -> None:
        self.assembly = assembly
        self.seqid = seqid
        self.source = source
        self.type = type
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.phase = phase
        self.attribute = attribute

    def to_dict(self):
        return {k: self.__getattribute__(k) for k in self.__slots__}


self = sg_object


# need to keep track of gene uids


for assembly_k, assembly_v in self.assemblies.items():
    for locus_k, locus_v in assembly_v.loci.items():
        for feature in locus_v.features:
            (assembly_k, locus_k, feature.start, feature.end, feature.strand)


a = Tabular(
    assembly=assembly_k,
    seqid=locus_k,
    source="sg",
    type=feature.type,
    start=feature.start,
    end=feature.end,
    strand=feature.strand,
)

inspect(a)

# an sg_object
# assembly order
# gene groups (not protein)


import json


class Clustermap(json.JSONEncoder):
    def __init__(self):
        self.uid = 0
        # holds mapping between sg objects and clustermap uids
        self.uid_dict = {}
    def _get_uid(self, obj=None):
        self.uid += 1
        self.uid_dict[str(self.uid)] = obj
        return str(self.uid)
    def _feature(self, feature_obj):
        return {
            "uid": self._get_uid(obj=feature_obj),
            "label": feature_obj.protein_id,
            "names": {
                "name": feature_obj.protein_id,
                "description": feature_obj.description,
            },
            "start": feature_obj.start,
            "end": feature_obj.end,
            "strand": feature_obj.strand,
        }
    def _locus(self, locus_name, locus_obj):
        return {
            "uid": self._get_uid(obj=locus_obj),
            "name": locus_name,
            "genes": [
                self._feature(feature_obj) for feature_obj in locus_obj.features
            ],
            "start": min([i.start for i in locus_obj.features]),
            "end": max([i.end for i in locus_obj.features]),
        }
    def _loci(self, assembly_obj):
        return [self._locus(locus_name=k, locus_obj=v) for k, v in assembly_obj.loci.items()]
    def _clusters(self, sg, assembly_order): 
        return {"clusters":[{"uid":self._get_uid(obj=sg.assemblies[k]), "name":k, "loci":self._loci(sg.assemblies[k])} for k in assembly_order]}
    def _get_gene_uid_of_protein_hash(self, protein_hash):
        return [k for k,v in self.uid_dict.items() if type(v).__name__ =="Feature" and v.protein_hash==protein_hash]
    def _flatten_list(self,x):
        return [item for row in x for item in row]
    def _create_group(self, k,v):
            return {"uid": self._get_uid(obj=k),"label":k, "genes":self._flatten_list([self._get_gene_uid_of_protein_hash(i) for i in v ])}
    def _create_groups_json(self, groupdict):
        return {"groups":[self._create_group(k,v) for k,v in groupdict.items()]}  
    def doit(self, sg, groupdict, assembly_order):
        return self._clusters(sg, assembly_order) | self._create_groups_json(groupdict)  | {"links":self._links(sg)}
    def _links(self, sg):
        sg_object.protein_comparison_to_df()
        res=list()
        # ugly because it first filters the protein_comparison table for proteins present in the sg_object
        for index, row in sg_object.protein_comparison.iterrows():
                for query_uid in self._get_gene_uid_of_protein_hash(row.query):
                    for target_uid in self._get_gene_uid_of_protein_hash(row.target):
                        res.append(
                            {"uid":self._get_uid(),
                                "query": {
                                    "uid": query_uid,
                                    "name": self._get_uid(),
                                },
                                "target": {
                                    "uid": target_uid,
                                    "name": self._get_uid(),
                                },
                                "identity": row.mod_score * 66,
                            }
                            # row["mod_score"] * 66 adjusts max mod score of 1.5 to ~100
                        )
        return res              


a=Clustermap()
#a.assembly_to_json(assembly_obj, bb)

z=a.doit(sg_object, groupdict=groupdict, assembly_order=list(sg_object.assemblies.keys()))


with open("/home/chase/Downloads/ttt/tempppp/clinker/clinker/plot/data.json", "w") as outfile:
    json.dump(z, outfile)
json.dumps(z)
BGC0000002

sg_object.annotate(use_neo4j_precalc=True)
sg_object.protein_comparison = []
sg_object.compare_proteins(append=True, cpus=20)
sg_object.protein_comparison_to_df()
sg_object.protein_comparison

assembly_order=["BGC0000002", "BGC0000001"]



def
