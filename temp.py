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


TRUTH_PROTEIN_UIDS = ["012774D5E4F15804", "02219101B8D12ECF", "18320CE1ED14CA4F", "19DF4D6715F3EF5A", "3AA6524AA5569D9A", "45302082D2D7430D", "49BAAA86ED17D693", "4A48029FA2F62DFC", "4F08A2D146BA8BBB", "8DABDB6FB9C444CE", "9030569BD80627CA", "91DA15E897A3C16C", "9A5CF4363385D6BD", "A22545A9337ACB7F", "A24055D4820B1374", "A50CF2CF1F40DEEC", "ACDA192F1F6F5B73", "AEB9CAA2344575D7", "B7B69593ECB95522", "D1A0B248EAF34BA3", "E231CD5F9EC69E7B", "E67BDB5B304AE16F", "FC06602FD6AEC428"]
# MATCH p=(n1 {external_id:"NC_015434.1"})-[r:ENCODES]->(p1)
# where r.start > 2562031 and r.end < 2617424
# RETURN collect(r.locus_tag)
TRUTH_LOCUS_TAGS =["VAB18032_RS11750", "VAB18032_RS11755", "VAB18032_RS11725", "VAB18032_RS11695", "VAB18032_RS11770", "VAB18032_RS31665", "VAB18032_RS11735", "VAB18032_RS11730", "VAB18032_RS11805", "VAB18032_RS11745", "VAB18032_RS11785", "VAB18032_RS11775", "VAB18032_RS11740", "VAB18032_RS11710", "VAB18032_RS11715", "VAB18032_RS11720", "VAB18032_RS11700", "VAB18032_RS11790", "VAB18032_RS11760", "VAB18032_RS11705", "VAB18032_RS11800", "VAB18032_RS11780", "VAB18032_RS11795"]


########################################
# Setup hmmsearch
########################################
h1=Path("/home/chase/Documents/socialgene_data/BGC0000001/socialgene_per_run/hmm_cache/socialgene_nr_hmms_file_with_cutoffs_1_of_1.hmm.gz")
h2=Path("/home/chase/Documents/socialgene_data/BGC0000001/socialgene_per_run/hmm_cache/socialgene_nr_hmms_file_without_cutoffs_1_of_1.hmm.gz")

h = HMMER()
h.hmmpress(h1)
h.hmmpress(h2)

########################################
# Read input BGC
########################################

sg_object = SocialGene()
sg_object.parse("/home/chase/Documents/github/kwan_lab/socialgene/sgpy/tests/python/data/bgc/BGC0000001.gbk.gz")
sg_object.annotate_proteins_with_neo4j(annotate_all=True)


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
# Calculate outdegree of HMMs
########################################

with GraphDriver() as db:
    _ = db.run("""
       MATCH (h1:hmm)-[r:ANNOTATES]->()
       WITH h1, count(r) as deg
       SET h1.outdegree = deg
     """)


########################################
# Prioritize proteins for search
########################################


log.info(f"{len(sg_object.proteins)} input proteins")

log.info(pd.DataFrame([{"protein":v.external_protein_id, "n_domains": len(v.domains)} for k,v in sg_object.proteins.items()]))
#log.info(pd.DataFrame([{"protein":k, "n_domains": len(v.domains)} for k,v in sg_object.proteins.items()]))


[{"protein":v.external_protein_id, "n_domains": len(v.domains)} for k,v in sg_object.proteins.items()]


so=np.argsort([len(v.domains) for k,v in sg_object.proteins.items()])


def sort_list_of_tuples(tup):
    arr = np.array(tup, dtype=[('col1', object), ('col2', int)])
    indices = np.argsort(arr['col2'])
    sorted_arr = arr[indices]
    sorted_tup = [(row['col1'], row['col2']) for row in sorted_arr]
    return sorted_tup

st = sort_list_of_tuples([(k, len(v.domains)) for k,v in sg_object.proteins.items()])


a=[]
doms=set()

for k,v in sg_object.proteins.items():
    for domain in v.domains:
        a.append({"protein_uid":k, "hmm_uid":domain.hmm_id})
        doms.add(domain.hmm_id)

a=pd.DataFrame(a)
a['n_hmms'] = a.groupby('protein_uid')['protein_uid'].transform('count')

no_dup = a.drop_duplicates()

with GraphDriver() as db:
    res = db.run("""
       WITH $doms as doms
       MATCH (h1:hmm)
       WHERE h1.uid in doms
       RETURN h1.uid as hmm_uid, h1.outdegree as outdegree
     """, doms=list(doms))
    z=res.to_df()

df=no_dup.merge(z, left_on="hmm_uid", right_on="hmm_uid")
df['total_outdegree'] = df.groupby('protein_uid')['outdegree'].transform('sum')
df =df.sort_values(by = ['total_outdegree'], ascending = [True], na_position = 'first')
df =df[['protein_uid', 'total_outdegree']].drop_duplicates()



def find_sim_protein(protein=sg_object.proteins['3AA6524AA5569D9A']):
    dv = list(set(protein.domain_vector))
    if dv:
        if len(dv) > 10:
            tol=2
        else:
            tol=1
        with GraphDriver() as db:
            res = db.run("""

        WITH $doms AS input_protein_domains
    WITH input_protein_domains, size(input_protein_domains) AS input_len
    MATCH (prot1:protein)<-[a1:ANNOTATES]-(h0:hmm)
    MATCH (prot1)<-[a2:ANNOTATES]-(h1:hmm)
    WHERE h0.uid IN input_protein_domains
    WITH input_len, prot1, count(DISTINCT (h0)) AS hmm_matches,  count(DISTINCT (h1)) AS all_annotations
    WHERE abs(input_len-hmm_matches) < $tol AND abs(input_len-all_annotations) < $tol
    MATCH (n1:nucleotide)-[e1:ENCODES]->(prot1)
    RETURN DISTINCT n1.uid as n, prot1.uid as p, input_len, hmm_matches, all_annotations, e1.start as start, e1.end as end
             """, doms=dv, tol=tol)
            return res.data()




df2 =pd.DataFrame({"in":[],"n":[], "p":[]})
bb=[]

for prot in list(df['protein_uid']):
    bb.extend([{"a":prot}|i for i in find_sim_protein(sg_object.proteins[prot])])

pd.DataFrame(bb)


