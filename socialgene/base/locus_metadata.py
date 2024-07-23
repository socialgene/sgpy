from collections import OrderedDict

from socialgene.utils.logging import log


class LocusAssemblyMetadata:
    # fmt: off
    __slots__= ["altitude","bio_material","bioproject","biosample","cell_line","cell_type","chromosome","clone","clone_lib","collected_by","collection_date","country","cultivar","culture_collection","db_xref","dev_stage","ecotype","environmental_sample","focus","germline","haplogroup","haplotype","host","identified_by","isolate","isolation_source","lab_host","lat_lon","macronuclear","map","mating_type","metagenome_source","mol_type","note","organelle","organism","pcr_primers","plasmid","pop_variant","proviral","rearranged","segment","serotype","serovar","sex","specimen_voucher","strain","sub_clone","submitter_seqid","sub_species","sub_strain","tissue_lib","tissue_type","transgenic","type_material","variety"]  # noqa: E231,E225
    # fmt: on

    def __init__(self, **kwargs) -> None:
        self.update(kwargs)

    def all_attributes(self):
        return OrderedDict({i: self.__getitem__(i) for i in self.__slots__})

    @property
    def values(self):
        return self.all_attributes().values()

    def __getitem__(self, item):
        try:
            return self.__getattribute__(item)
        except Exception as e:
            log.debug(e)
            return None

    # add a function that takes a dict as input and sets any attributes in the dict keys
    def update(self, d: dict):
        for k, v in d.items():
            temp_v = v
            if isinstance(v, list):
                # if info is a list, collapse into a single string
                if len(v) > 1:
                    temp_v = ";".join(v)
                else:
                    temp_v = v[0]
            if not isinstance(temp_v, str):
                temp_v = str(temp_v)
            temp_v.replace("\t", "")
            temp_v.replace("\n", "")
            temp_v.replace("\r", "")
            temp_v.strip()
            if k in self.__slots__:
                setattr(self, k, temp_v)
