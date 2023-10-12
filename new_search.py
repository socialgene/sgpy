import time
from socialgene.search.hmmer import SearchDomains
from socialgene.base.socialgene import SocialGene

now = time.time()


input_gbk_path = "/home/chase/Documents/data/mibig/3_1/mibig_gbk_3.1/BGC0001850.gbk"
hmm_dir = None
hmm_dir = (
    "/home/chase/Documents/socialgene_data/streptomyces/socialgene_per_run/hmm_cache"
)


a = SocialGene()
a.parse(input_gbk_path)

search_object = SearchDomains(
    gbk_path=input_gbk_path,
    hmm_dir=hmm_dir,
    assemblies_must_have_x_matches=10,
    nucleotide_sequences_must_have_x_matches=10,
    gene_clusters_must_have_x_matches=10,
)


search_object.outdegree_table()

search_object.prioritize_input_proteins(
    max_domains_per_protein=1,
    max_outdegree=None,
    max_query_proteins=20,
)
search_object.outdegree_table()


search_object.search_db()

search_object.process_results()


search_object.write_clustermap(
    "/home/chase/Downloads/clinker-master(2)/clinker-master/clinker/plot/data.json"
)

later = time.time()
int(later - now)
