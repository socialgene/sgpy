import time

from socialgene.clustermap.serialize import SerializeToClustermap
from socialgene.search.hmmer import SearchDomains

start_time = time.time()
# input_gbk_path = "/home/chase/Downloads/para/a/GCF_002362315.1_ASM236231v1_genomic.gbff.gz"
input_gbk_path = "/home/chase/Documents/data/mibig/3_1/mibig_gbk_3.1/BGC0001646.gbk"
hmm_dir = None
# hmm_dir = "/home/chase/Downloads/meh/socialgene_per_run/hmm_cache"
hmm_dir = "/home/chase/Downloads/para/socialgene_per_run/hmm_cache"


search_object = SearchDomains(
    gbk_path=input_gbk_path,
    hmm_dir=hmm_dir,
    use_neo4j_precalc=True,
    assemblies_must_have_x_matches=0.6,
    nucleotide_sequences_must_have_x_matches=0.6,
    gene_clusters_must_have_x_matches=0.6,
    break_bgc_on_gap_of=10000,
    target_bgc_padding=15000,
)


# search_object.outdegree_table


search_object.prioritize_input_proteins(
    max_domains_per_protein=2,
    max_outdegree=100000,
    max_query_proteins=5,
    scatter=False,
    # bypass_locus=["MicB006_2899"]
    # bypass_pid=["AXA20096.1"],
)

# search_object.outdegree_table

# TODO frac is redundant with the outdegree table
search_object.search(only_culture_collection=False, frac=0.75)


search_object.filter()
search_object.label_clusters()


search_object._bgc_regions_to_sg_object(search_object._primary_bgc_regions())

_ = search_object.sg_object.annotate_proteins_with_neo4j(
    protein_uids=None, annotate_all=True, progress=False
)


df = search_object._primary_bgc_regions()


for i, row in df.iterrows():
    search_object.sg_object.assemblies[row["assembly_uid"]].get_locus_by_uid(
        row["nucleotide_uid"]
    ).add_bgcs_by_start_end(
        start=row["n_start"],
        end=row["n_end"],
        uid=row["cluster"],
    )


temp = search_object.input_assembly.loci[search_object.input_bgc_id]
temp.add_bgcs_by_feature(features=temp.features)

search_object._create_links()
search_object._choose_group()


z = SerializeToClustermap(
    sg_object=search_object.sg_object,
    bgc_order=list(search_object.sg_object.assemblies.keys()),
    link_df=search_object.link_df,
    group_df=search_object.group_df,
)

z.write("/home/chase/Downloads/clinker-master(2)/clinker-master/clinker/plot/data.json")
print("--- %s seconds ---" % (time.time() - start_time))

# search_object.rich_table(search_object.user_friendly_hit_df())


# self = search_object
