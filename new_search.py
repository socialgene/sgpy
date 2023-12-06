# The code is performing a series of operations using the SocialGene package and the SearchDomains
# class from the socialgene.search.hmmer module. Here is a breakdown of what the code is doing:

from socialgene.base.socialgene import SocialGene
from socialgene.clustermap.serialize import SerializeToClustermap
from socialgene.search.hmmer import SearchDomains
from socialgene.config import env_vars

env_vars["NEO4J_URI"] = "bolt://localhost:7687"

input_gbk_path = "/home/chase/Documents/data/mibig/3_1/mibig_gbk_3.1/BGC0001850.gbk"
# input_gbk_path = "/home/chase/Documents/data/mibig/3_1/mibig_gbk_3.1/BGC0001848.gbk"
hmm_dir = None
# hmm_dir = "/home/chase/Downloads/meh/socialgene_per_run/hmm_cache"
hmm_dir = "/home/chase/Documents/socialgene_data/v0_4_1/streptomyces/socialgene_per_run/hmm_cache"

a = SocialGene()
a.parse(input_gbk_path)
len(a.proteins)

search_object = SearchDomains(
    gbk_path=input_gbk_path,
    hmm_dir=hmm_dir,
    use_neo4j_precalc=True,
    assemblies_must_have_x_matches=0.6,
    nucleotide_sequences_must_have_x_matches=0.6,
    gene_clusters_must_have_x_matches=0.6,
    break_bgc_on_gap_of=10000,
    target_bgc_padding=50000,
)
self = search_object

search_object.outdegree_table


search_object.prioritize_input_proteins(
    max_domains_per_protein=5,
    max_outdegree=None,
    max_query_proteins=8,
    scatter=False,
    # loci=["MicB006_2899"]
    # bypass_eid=["AXA20096.1"],
)

search_object.outdegree_table

# TODO frac is redundant with the outdegree table
search_object.search(only_culture_collection=True, frac=0.75)


search_object.filter()
search_object.label_clusters()


search_object._bgc_regions_to_sg_object(self._primary_bgc_regions())

_ = self.sg_object.annotate_proteins_with_neo4j(
    protein_uids=None, annotate_all=True, progress=False
)


df = self._primary_bgc_regions()


for i, row in df.iterrows():
    self.sg_object.assemblies[row["assembly_uid"]].get_locus_by_uid(
        row["nucleotide_uid"]
    ).add_bgcs_by_start_end(
        start=row["n_start"],
        end=row["n_end"],
        uid=row["cluster"],
    )


temp = self.input_assembly.loci[self.input_bgc_id]

temp.add_bgcs_by_feature(features=temp.features)


self._create_links()
self._choose_group()


z = SerializeToClustermap(
    sg_object=self.sg_object,
    bgc_order=list(self.sg_object.assemblies.keys()),
    link_df=self.link_df,
    group_df=self.group_df,
)

z.write("/home/chase/Downloads/clinker/clinker/plot/data.json")


# for v in self.sg_object.proteins.values():
#     v.add_sequences_from_neo4j()

# self.sg_object.assemblies["GCF_000720145.1"].loci["NZ_JOIS01000040.1"].write_genbank(
#     "test.gbk", start=34872, end=54001
# )

# inspect(self.sg_object.assemblies["GCF_000720145.1"])
