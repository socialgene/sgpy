import json

from socialgene.utils.logging import log


class SerializeToClustermap:
    """
    Take a sg_object and serialize all BGC object to clustermap.js format
    """

    def __init__(self, sg_object, bgc_order, link_df, group_df):
        self._sg_object = sg_object
        self._bgc_order = bgc_order
        self._input_assembly = self._sg_object.assemblies[bgc_order[0]]
        self._link_df = link_df
        self._group_df = group_df
        self._reset()

    def _reset(self):
        self._uid = 0
        self._uid_dict = {}
        self.feature_to_cmap_uid_dict = {}  # {"feature_obj": "clustermap_uid"}

    def _flatten_list(x):
        return sorted([item for row in x for item in row], reverse=True)

    def _get_uid(self, obj=None):
        # increment then return a uid as a string
        self._uid += 1
        self._uid_dict[str(self._uid)] = obj
        return str(self._uid)

    def _build(self):
        self._reset()
        return self._clusters() | self._create_links_dict() | self._create_groups_dict()

    def _clusters(self):
        log.info("Creating clustermap.js clusters")
        return {
            "clusters": [
                {
                    "uid": self._get_uid(obj=self._sg_object.assemblies[k]),
                    "name": self._sg_object.assemblies[k].name,
                    "loci": self._loci(self._sg_object.assemblies[k]),
                }
                for k in self._bgc_order
            ]
        }

    def _feature(self, feature_obj):
        feat_id = self._get_uid(obj=feature_obj)
        self.feature_to_cmap_uid_dict.update({feature_obj: feat_id})
        return {
            "uid": feat_id,
            "label": feature_obj.external_id,
            "names": {
                "name": feature_obj.external_id,
                "description": feature_obj.description,
                "id": feature_obj.external_id,
                "locus_tag": feature_obj.locus_tag,
            },
            "start": feature_obj.start,
            "end": feature_obj.end,
            "strand": feature_obj.strand,
        }

    def _locus(self, locus_name, locus_obj):
        return {
            "uid": self._get_uid(obj=locus_obj),
            "name": locus_name,
            "genes": [self._feature(feature_obj) for feature_obj in locus_obj.features],
            "start": min((i.start for i in locus_obj.features)),
            "end": max((i.end for i in locus_obj.features)),
        }

    def _loci(self, assembly_obj):
        # return [
        #     self._locus(locus_name=k, locus_obj=v) for k, v in assembly_obj.loci.items()
        # ]
        return [
            self._locus(locus_name=k, locus_obj=gc)
            for k, v in assembly_obj.loci.items()
            for gc in v.gene_clusters
        ]

    def _create_groups_dict(self):
        """
        Create a dictionary of clustermap groups.

        Returns:
            A dictionary containing a list of groups, where each group is represented as a dictionary
            with keys "uid", "label", and "genes". e.g. {"groups": [{"uid": "1", "label": "group1", "genes": ["1", "2", "3"]}]
        """
        # Look in the clustermap uid dict (feature_to_cmap_uid_dict) for the cmap uid of each feature (query & target) in each group
        self._group_df["query_cmap_uid"] = self._group_df["query"].apply(
            lambda x: self.feature_to_cmap_uid_dict.get(x)
        )
        self._group_df["target_cmap_uid"] = self._group_df["target"].apply(
            lambda x: self.feature_to_cmap_uid_dict.get(x)
        )
        groups = (
            self._group_df.groupby(["target_cmap_uid", "target"])["query_cmap_uid"]
            .apply(list)
            .reset_index()
        )
        return {
            "groups": [
                {
                    "uid": self._get_uid(
                        obj=f"{i.target.external_id} {i.target.description}"
                    ),
                    "label": f"{i.target.external_id} {i.target.description}",
                    "genes": i.query_cmap_uid,
                }
                for x, i in groups.iterrows()
            ]
        }

    def _create_links_dict(self):
        log.info("Creating clustermap.js links")

        return {
            "links": [
                {
                    "uid": self._get_uid(obj=None),
                    "target": {
                        "uid": self.feature_to_cmap_uid_dict[i["target"]],
                        "name": self.feature_to_cmap_uid_dict[i["target"]],
                    },
                    "query": {
                        "uid": self.feature_to_cmap_uid_dict[i["query"]],
                        "name": self.feature_to_cmap_uid_dict[i["query"]],
                    },
                    "identity": i.score,
                }
                for x, i in self._link_df.iterrows()
                if i["query"] in self.feature_to_cmap_uid_dict
                and i["target"] in self.feature_to_cmap_uid_dict
            ]
        }

    def write(self, outpath):
        log.info(f"Writing clustermap.js output to: {outpath}")
        with open(outpath, "w") as outfile:
            json.dump(
                self._build(),
                outfile,
            )
