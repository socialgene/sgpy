# MATCH p=(n1 {external_id:"NC_015434.1"})-[r:ENCODES]->(p1)
# where r.start > 2562031 and r.end < 2617424
# RETURN collect(p1.uid)
from collections import defaultdict
from socialgene.hmm.hmmer import HMMER
from socialgene.base.socialgene import SocialGene
from rich import inspect
from pathlib import Path
from socialgene.neo4j.neo4j import GraphDriver
from socialgene.utils.logging import log
import pandas as pd
import numpy as np
from copy import deepcopy
from datetime import datetime
import logging

logging.getLogger("neo4j").setLevel(logging.WARNING)
logging.getLogger().setLevel(logging.INFO)
# fmt:off

# Do this in BASH first then
#       curl https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/408/515/GCF_003408515.1_ASM340851v1/GCF_003408515.1_ASM340851v1_genomic.gbff.gz > /tmp/GCF_003408515.1_ASM340851v1_genomic.gbff.gz
#
#       outdir='/home/chase/Documents/socialgene_data/BGC0001848'
#       outdir_download_cache='/tmp/socialgene_data/cache'
#
#       nextflow run socialgene/sgnf -r 502e371b77 \
#           -profile docker \
#           --outdir $outdir \
#           --local_genbank '/tmp/GCF_003408515.1_ASM340851v1_genomic.gbff.gz' \
#           --outdir_download_cache $outdir_download_cache \
#           --mmseqs_steps '90 70 50' \
#           --blastp \
#           --hmmlist 'antismash' \
#           --build_database
#

# MATCH p=(n1 {external_id:"NZ_CP030865.1"})-[r:ENCODES]->(p1)
# where r.start >= 3038865 and r.end <= 3097062
# RETURN collect(p1.uid)
TRUTH_PROTEIN_UIDS = ["zwG8iK66IVseCDcI3OwBT-tNImw4Hbu7", "zX5ILAQELaLV0pDShayy7nx70loP2LLZ", "z47L7LvD5PZmYwKXEOn6Q5-ksyuSDXSY", "xT0WaclZ6XPn9WwTNrO5yYVy3IyXEE4E", "w8z5aogfz7vLz6TiPt2Wt2pWxSCfl6_C", "uNrCjGQbe9LsL9umHQWqvG2OQykvFxFw", "uMvE-xgz2R-SI_mzY7HucLRXLm2WkwYC", "u6iVhXiaP7LY-ncGV0nMtvrJlZD7vjsx", "qEcFCOXPHpf_D6FiGPea8DkV_AWWZjMy", "orV-86yt1TbEwWGDAkU3bVNpe459Hl1B", "oS3H0Kb3lOQz6-wnh14vTmcPhT2_FDT_", "lsxvnXk8kxp2OelIap6TcoguSBVOEMnk", "kZ3wvfvDYav1WlgZPVfhYrfPQsuIHUmn", "jb7i54nlRxscCRkVc-T3G7DhAvzwdPYW", "etULepIwK3A8cTlCpg7R-CAWJISYYy0t", "ennTHgEojAJicxoKyyb-1MWJz5q0YNGv", "aX40NgOFtJ7diWuBmRFVm-6Lm16iECFu", "_bM-9OJARu4IC1V8jchTXDIaTPuLUhCe", "_2LS290Ljg7VSPbsUSwDel51WcQzbxNM", "_0OhI8CwFqgCa97qVet-jbUA7WN2Q3XF", "X5QhwaU2vttuwby1ZKaJj8jAkDySa8SS", "UpM8jOJT2xvzuU9Nbq6foihOdCQg0JR0", "U3Cgn6zYdFSQVEZFTRUnAo1zaJfp9-Q4", "TPii9Xdnw4CEdcZiTAt7qkdVP0ZQZAG8", "TJEr78TMOE8Dt84sRoNjM4Gr05mHQPaq", "SZWGdb4CRIEdkRg6gcygIGi5mubgmLjz", "Rie2hnRZGyQbrC0oV6vKmK4xyx8UM3c6", "QJgl69uU344A9frRrwPiY9jOgFPJehcg", "O6l97fsTY-BTUyd3eZ0kcxotvrBKJeeZ", "M-XQN-LZbGTFgj0rp2dFnIsF9bQm9chC", "M-VPo39ep4Utq6_Fru4JuNM9rp3Wmu45", "KbJwxUWA9hBamAnz7boyic9xuL7YpdFp", "IU-V8ZvWVd1C1Z-evwLfP_KkJT7mZCvA", "Ft0znrFotJAJw4O3UQcsSde2fg7joffB", "FPKQ-PvVfRN5qj31YHB_fDT2Q71Y1bLw", "Dtbpljoi8PkqXYUIh1uwfB3uyGtIvIYz", "D8hihmbo9jMJDf3zV2YsgAKSTgxHwc05", "C4J0nKGgGspm3rM_DkeKDTZqvqQSCKjS", "BIkruFDUoKjckgk0j6oZUrBk_oAj16YG", "9FIPGbtPdODtNGOjE2HrlxtNpnvB7dOS", "6cyDFNPg9PUuv6MoerdrWu2fEKfbz-dP", "6DR3nKfuJiYk9MzbWL6Ly_m0HrRkiTBo", "4tpU-EDbIbPxLJPJeHH3FT034Y2dK1it", "4EeqbFY1V20e1Z8RNjccxJsJaspnPGB1", "45J3u_EljbjWPJr9yQio43TzW0K6Wa3O", "0yVaZrn7OtVTYy5bfiHLRUYCRvmhh5QF", "0suO2ImDjlLb-G5pYR1evP6CDaM4Kmma", "0R8VEspGEA3N9FuEVEymKuXmWrmhIbEc", "-gx0Evqh419j7dgxxuecybcdfqoH5Zjp", "-Kk5iCbSLfpI6pXAOskPjSrFhFL4iGBH"]

# MATCH p=(n1 {external_id:"NZ_CP030865.1"})-[r:ENCODES]->(p1)
# where r.start => 3038865 and r.end =< 3097062
# RETURN collect(r.locus_tag)
TRUTH_LOCUS_TAGS =["MicB006_RS14515", "MicB006_RS14450", "MicB006_RS14455", "MicB006_RS14495", "MicB006_RS14490", "MicB006_RS14560", "MicB006_RS14430", "MicB006_RS14480", "MicB006_RS14375", "MicB006_RS14600", "MicB006_RS14420", "MicB006_RS14615", "MicB006_RS14385", "MicB006_RS14505", "MicB006_RS14550", "MicB006_RS14580", "MicB006_RS14570", "MicB006_RS14545", "MicB006_RS14485", "MicB006_RS14400", "MicB006_RS14475", "MicB006_RS14405", "MicB006_RS14610", "MicB006_RS14470", "MicB006_RS14595", "MicB006_RS14535", "MicB006_RS14520", "MicB006_RS14415", "MicB006_RS14435", "MicB006_RS14620", "MicB006_RS14565", "MicB006_RS14410", "MicB006_RS14590", "MicB006_RS14500", "MicB006_RS14445", "MicB006_RS14460", "MicB006_RS14555", "MicB006_RS14525", "MicB006_RS14585", "MicB006_RS14540", "MicB006_RS14530", "MicB006_RS14465", "MicB006_RS14440", "MicB006_RS14605", "MicB006_RS14575", "MicB006_RS14425", "MicB006_RS14380", "MicB006_RS14510", "MicB006_RS14390", "MicB006_RS14395"]

then  = datetime.now()

########################################
# Setup hmmsearch
########################################
h1=Path("/home/chase/Documents/socialgene_data/BGC0001848/socialgene_per_run/hmm_cache/socialgene_nr_hmms_file_with_cutoffs_1_of_1.hmm.gz")
h2=Path("/home/chase/Documents/socialgene_data/BGC0001848/socialgene_per_run/hmm_cache/socialgene_nr_hmms_file_without_cutoffs_1_of_1.hmm.gz")

h = HMMER()
h.hmmpress(h1)
h.hmmpress(h2)

########################################
# Read input BGC
########################################

sg_object = SocialGene()
#sg_object.parse("/tmp/BGC0001848.gbk")
sg_object.parse("/home/chase/Documents/data/mibig/3_1/mibig_gbk_3.1/BGC0001646.gbk")
_ = sg_object.annotate_proteins_with_neo4j(annotate_all=True)



########################################
# Annotate input BGC with the HMMs
########################################

# h1.with_suffix("") removes the ".gz" extension
sg_object.annotate_with_hmmscan(
    protein_id_list=list(sg_object.proteins.keys()), hmm_filepath=h1.with_suffix(""), cpus=1, cut_ga=True
)
sg_object.annotate_with_hmmscan(
    protein_id_list=list(sg_object.proteins.keys()), hmm_filepath=h2.with_suffix(""), cpus=1
)

########################################
# Calculate outdegree of HMMs and proteins
########################################

# Commented out because it modifies the DB, so only need to run this once

# with GraphDriver() as db:
#     _ = db.run("""
#        MATCH (h1:hmm)-[r:ANNOTATES]->()
#        WITH h1, count(r) as deg
#        SET h1.outdegree = deg
#      """)

# with GraphDriver() as db:
#     _ = db.run("""
#        MATCH (p1:protein)<-[:ENCODES]-(n1:nucleotide)
#        WITH p1, count(DISTINCT n1) as deg
#        SET p1.nucleotide_outdegree = deg
#      """)

########################################
# Prioritize proteins for search
########################################

# log.info(f"{len(sg_object.proteins)} input proteins")
# log.info(pd.DataFrame([{"protein":v.external_protein_id, "n_domains": len(v.domains)} for k,v in sg_object.proteins.items()]))
# log.info(pd.DataFrame([{"protein":k, "n_domains": len(v.domains)} for k,v in sg_object.proteins.items()]))
# [{"protein":v.external_protein_id, "n_domains": len(v.domains)} for k,v in sg_object.proteins.items()]

########################################
# Prioritize input proteins for search
########################################

# Prioritize based in the cumulative outdegree of HMM nodes per input protein
input_protein_domain_df=[]
doms=set()

for k,v in sg_object.proteins.items():
    for domain in v.domains:
        input_protein_domain_df.append({"protein_uid":k, "hmm_uid":domain.hmm_id})
        doms.add(domain.hmm_id)

input_protein_domain_df=pd.DataFrame(input_protein_domain_df)
input_protein_domain_df['n_hmms'] = input_protein_domain_df.groupby('protein_uid')['protein_uid'].transform('count')

input_protein_domain_df = input_protein_domain_df.drop_duplicates()

with GraphDriver() as db:
    res = db.run("""
       WITH $doms as doms
       MATCH (h1:hmm)
       WHERE h1.uid in doms
       RETURN h1.uid as hmm_uid, h1.outdegree as outdegree
     """, doms=list(doms))
    z=res.to_df()

input_protein_domain_df=input_protein_domain_df.merge(z, left_on="hmm_uid", right_on="hmm_uid")
input_protein_domain_df['total_outdegree'] = input_protein_domain_df.groupby('protein_uid')['outdegree'].transform('sum')
input_protein_domain_df =input_protein_domain_df.sort_values(by = ['total_outdegree'], ascending = [True], na_position = 'first')
input_protein_domain_df =input_protein_domain_df[['protein_uid', 'total_outdegree']].drop_duplicates()

# Take the top 50% of input proteins that have the lowest outdegree
#input_protein_domain_df =  input_protein_domain_df.head(int(len(input_protein_domain_df)/3))
input_protein_domain_df=input_protein_domain_df.head(5)
########################################
# Find DB proteins that match the input proteins
########################################

def find_sim_protein(protein):
    # This could be done for the whole BGC at once, but is done protein by protein to be more atomic
    dv = list(set(protein.domain_vector))[0:3]
    with GraphDriver() as db:
            res = db.run("""
    WITH $doms AS input_protein_domains, size($doms) AS input_len
    MATCH (prot1:protein)<-[a1:ANNOTATES]-(h0:hmm)
    WHERE h0.uid IN input_protein_domains
    MATCH (prot1)<-[a2:ANNOTATES]-(h1:hmm)
    WITH input_len, prot1, count(DISTINCT (h0)) AS hmm_matches,  count(DISTINCT (h1)) AS all_annotations
    WHERE abs(input_len-hmm_matches)  AND abs(input_len-all_annotations) < $tol
    MATCH (n1:nucleotide)-[e1:ENCODES]->(prot1)
    RETURN DISTINCT n1.uid as n, prot1.uid as p, input_len, hmm_matches, all_annotations, e1.start as start, e1.end as end
             """, doms=dv, tol=tol)
            return res.data()




bb=[]

for prot in list(input_protein_domain_df['protein_uid']):
    bb.extend([{"a":prot}|i for i in find_sim_protein(sg_object.proteins[prot])])

bb=pd.DataFrame(bb)

zz =bb.sort_values(by = ['start'], ascending = [True], na_position = 'first')
zz.reset_index(inplace=True, drop=True)

start_series = deepcopy(zz.start)


# groups = []
# for index, row in start_series.iterrows():
#     raise
#     index_plus_one = index + 1
#     in_window = (zz.loc[index:, 'start'] - row.start) < 10000
#     groups.append(zz.index[in_window])


# np.diff(zz.start)

# ind = np.where(np.diff(zz.start))[0]

# np.split(start_series, np.where(np.diff(start_series) != 10000)[0]+1)


########################################
# Number of input proteins that had a hit in a nucleotide sequence
########################################
prioritized_nucleotide_seqs = zz.filter(['n','a']).drop_duplicates().groupby(['n']).count()
prioritized_nucleotide_seqs = prioritized_nucleotide_seqs.sort_values(by = ['a'], ascending = [False], na_position = 'first')
prioritized_nucleotide_seqs.reset_index(inplace=True, drop=False)
prioritized_nucleotide_seqs.rename(columns={"n":"nucleotide_uid", "a": "count_of_matched_inputs"}, inplace=True)

prioritized_nucleotide_seqs=prioritized_nucleotide_seqs[prioritized_nucleotide_seqs["count_of_matched_inputs"] > len(input_protein_domain_df) * 0.75]

prioritized_nucleotide_seqs

#################




def find_best_cluster(df, tolerance= 40000):
    # Within an input input_protein_domain_df, find groups of genes within a tolerance of the number of
    # AAs from each other's starts
    # and return a input_protein_domain_df representing the group with the hits to the largest number of distinct input proteins
    s=0
    e=0
    max_match = (0,0,0)
    for i,diff in enumerate(np.diff(df.start)):
        if diff > tolerance:
            # len(set(df.iloc[s:e].a)))  count how many input proteins matched within a window
            if (matches := len(set(df.iloc[s:e].a))) > max_match[2]:
                max_match = (s, e, matches)
            # start new window
            s = i + 1
            e = s
        else:
            # expand window
            e = i
    # account for no splits
    if max_match ==  (0,0,0):
        max_match = (s, e, len(set(df.iloc[s:e].a)))
    return df.iloc[max_match[0]:max_match[1]]

bro=[]
for name, group in zz[zz['n'].isin(list(prioritized_nucleotide_seqs.nucleotide_uid))].groupby(['n']):
    bro.append(find_best_cluster(group))

pd.concat(bro)

# BGC0001848      G3Abi4t2KLNuJ2MzbteeLbGZ3-sAQSR3
# B006            kNauYNQ_TytXuMSFT8jxDuocmtJzCNQk


"G3Abi4t2KLNuJ2MzbteeLbGZ3-sAQSR3" in set(pd.concat(bro).n)
"kNauYNQ_TytXuMSFT8jxDuocmtJzCNQk" in set(pd.concat(bro).n)
"5u1SlwuuHRJrhUc12vHsG4E2CUPcrQZa" in set(pd.concat(bro).n)


now  = datetime.now()
duration = now - then
duration_in_s = duration.total_seconds()

minutes = divmod(duration_in_s, 60)[0]
