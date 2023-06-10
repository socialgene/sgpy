import itertools
import json
import uuid

from socialgene.base.compare_protein import CompareProtein
from socialgene.utils.np_json_converter import np_json_converter

# This code is meant to be used to turn a SocialGene object into json
# that can be read by https://github.com/gamcil/clustermap.js

# step 1: create uuids for everything that will go into clustermap
# step 2: create clusters

# done for any comparisons
# step 3: create links
# step 4: crete groups

# should be generalized so other classes/methods can:
# create a clustermap for a single socialgene object
# create a clustermap for multiple socialgene objects
# create a clustermap for comparative socialgene objects
# add to existing clustermaps


# TODO: error, uid for protein must be assigned within locus because duplicates can exist in different/same loci


class UuidCount:
    """Holds UUIDs for clustermap.js object"""

    __slots__ = [
        "uuid_counter",
    ]

    def __init__(
        self,
    ):
        self.uuid_counter = 0


class ClustermapUuids(UuidCount):
    """For flexibility `sg_object` was purposefully provided as an argument to each function and not as a class variable

    Returns:
        _type_: _description_
    """

    # track = UuidCount()  # track uuids across all instances of ClustermapUuids()

    def __init__(
        self,
    ):
        super(ClustermapUuids, self).__init__()
        self.uuid_dict = {}

    @staticmethod
    def build_protein_key(assembly_key, locus_key, locus_value):
        return f"{assembly_key}_{locus_key}_{locus_value}"

    def create_clustermap_uuids(self, sg_object):
        # create protein ids per-assembly/locus to handle nr proteins (clustermap ids must be unique even for nr)
        for assembly_key, assembly_values in sg_object.assemblies.items():
            # create uuid for assebmbly id
            if assembly_key not in self.uuid_dict:
                self.uuid_dict[assembly_key] = f"uuid_{str(self.uuid_counter)}"
                self.uuid_counter += 1
            for loci_key, loci_value in assembly_values.loci.items():
                # create uuid for locus id
                if loci_key not in self.uuid_dict:
                    self.uuid_dict[loci_key] = f"uuid_{str(self.uuid_counter)}"
                    self.uuid_counter += 1
                # create uuid for protein (id)
                for feature in loci_value.features:
                    protein_key = self.build_protein_key(
                        assembly_key=assembly_key,
                        locus_key=loci_key,
                        locus_value=feature.protein_hash,
                    )
                    if protein_key not in self.uuid_dict:
                        self.uuid_dict[protein_key] = f"uuid_{str(self.uuid_counter)}"
                        self.uuid_counter += 1


class Clustermap(ClustermapUuids, CompareProtein):
    def __init__(self, primary_assembly=None):
        super().__init__()
        self.primary_assembly = primary_assembly
        self.clusters = []
        self.groups = []
        self.links = []
        self.prot_loc = {}

    def write_clustermap(self, outpath: str = None):
        with open(outpath, "w") as handle:
            handle.write(self.dumps())

    def dumps(self):
        # np_json_converter() to fix: TypeError: Object of type int64 is not JSON serializable
        return json.dumps(
            {
                "clusters": self.clusters,
                "links": self.links,
                "groups": self.groups,
            },
            default=np_json_converter,
        )

    def add_cluster(self, sg_object):
        for assembly_k, assembly_v in sg_object.assemblies.items():
            _loci = []
            for locus_k, locus_v in assembly_v.loci.items():
                _genes = []
                _gene_starts = []
                _gene_ends = []
                for feature in locus_v.features:
                    if feature.feature_is_protein():
                        protein_key = self.build_protein_key(
                            assembly_key=assembly_k,
                            locus_key=locus_k,
                            locus_value=feature.protein_hash,
                        )
                        # Keep track of {protein: protein_key}
                        if feature.protein_hash not in self.prot_loc:
                            self.prot_loc[feature.protein_hash] = []
                        self.prot_loc[feature.protein_hash].append(protein_key)
                        _genes.append(
                            {
                                "uid": self.uuid_dict.get(protein_key),
                                "label": sg_object.proteins[
                                    feature.protein_hash
                                ].external_protein_id,
                                "names": {
                                    "name": sg_object.proteins[
                                        feature.protein_hash
                                    ].external_protein_id,
                                    "description": sg_object.proteins[
                                        feature.protein_hash
                                    ].description,
                                },
                                "start": feature.start,
                                "end": feature.end,
                                "strand": feature.strand,
                            }
                        )
                        _gene_starts.append(feature.start)
                        _gene_ends.append(feature.end)
                if not _gene_starts:
                    # stop here if no genes
                    break
                _loci.append(
                    {
                        "uid": self.uuid_dict.get(locus_k),
                        "name": locus_k,
                        "start": min(_gene_starts),
                        "end": max(_gene_ends),
                        "genes": _genes,
                    }
                )
            self.clusters.append(
                {
                    "uid": self.uuid_dict.get(assembly_k),
                    "name": assembly_k,
                    "loci": _loci,
                }
            )

    def add_groups(self, sg_object, cutoff: int = 0):
        # pandas explanaton- group pandas df and return as dict = {query:[target, target, target]}
        sg_object.protein_comparison_to_df()
        for k, v in (
            sg_object.protein_comparison.loc[
                sg_object.protein_comparison["mod_score"] > cutoff
            ]
            .groupby("query", group_keys=False)["target"]
            .apply(list)
            .to_dict()
            .items()
        ):
            # one list containing query and target protein hash_ids
            v.append(k)
            # add query protein
            gene_list = set()
            # see below for explanation of this list comprehension
            temp = [val for key, val in self.uuid_dict.items() if key.endswith(k)]
            for i in temp:
                gene_list.add(i)
            del temp
            # add matched proteins
            for i in set(v):
                # to deal with nr proteins having unique ids for clustermap,
                # the proteins uid keys are concatenated ids for assembly, locus, protein
                # so we have to find all uuid_dict keys that end with the protein hash
                # then those uuids to the list
                for ii in [
                    val for key, val in self.uuid_dict.items() if key.endswith(i)
                ]:
                    gene_list.add(ii)
            # label is the query protein
            # genes are match proteins + query protein
            if k in sg_object.proteins:
                if (
                    sg_object.proteins[k].external_protein_id is not None
                    and sg_object.proteins[k].description is not None
                ):
                    label = f"{sg_object.proteins[k].external_protein_id}__{sg_object.proteins[k].description}"
                elif sg_object.proteins[k].external_protein_id is not None:
                    label = sg_object.proteins[k].external_protein_id
                elif sg_object.proteins[k].description is not None:
                    label = sg_object.proteins[k].description
                else:
                    label = "None"
                self.groups.append(
                    {
                        "uid": str(uuid.uuid4()),
                        "label": label,
                        "genes": list(gene_list),
                        "hidden": False,
                    }
                )

    def add_links(self, sg_object):
        sg_object.protein_comparison_to_df()
        # ugly because it first filters the protein_comparison table for proteins present in the sg_object
        for index, row in sg_object.protein_comparison.iterrows():
            if row[0] in self.prot_loc and row[1] in self.prot_loc:
                for x, y in itertools.product(
                    self.prot_loc[row[0]],
                    self.prot_loc[row[1]],
                ):
                    self.links.append(
                        {
                            "uid": str(uuid.uuid4()),
                            "query": {
                                "uid": self.uuid_dict.get(x),
                                "name": row[0],
                            },
                            "target": {
                                "uid": self.uuid_dict.get(y),
                                "name": row[1],
                            },
                            "identity": row["mod_score"] * 66,
                        }
                        # row["mod_score"] * 66 adjusts max mod score of 1.5 to ~100
                    )
