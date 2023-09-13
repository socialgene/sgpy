import logging
from pathlib import Path

from socialgene.base.socialgene import SocialGene
from socialgene.compare_proteins.hmmer import CompareDomains
from socialgene.hmm.hmmer import HMMER

logging.getLogger("neo4j").setLevel(logging.WARNING)
logging.getLogger().setLevel(logging.INFO)


gbk_path = "/home/chase/Documents/data/mibig/3_1/mibig_gbk_3.1/BGC0001850.gbk"
hmm_dir = "/home/chase/Downloads/old3/work/06/45cc132a2e746452c578fefc8c34a7"
# hmm file with cutoffs
h1 = Path(hmm_dir, "socialgene_nr_hmms_file_with_cutoffs_1_of_1.hmm.gz")
# hmm file without cutoffs
h2 = Path(hmm_dir, "socialgene_nr_hmms_file_without_cutoffs_1_of_1.hmm.gz")
########################################
# Setup hmmsearch
########################################
# hmmpress
HMMER().hmmpress(h1)
HMMER().hmmpress(h2)

########################################
# Read input BGC
########################################

sg_object = SocialGene()

# parse input BGC
sg_object.parse(gbk_path)

# fragile implementation; add string input BGC id to keep track/separate input BGC from those added to class later
input_bgc_name = list(sg_object.assemblies.keys())[0]
modified_input_bgc_name = f"socialgene_query_{input_bgc_name}"

# change input assembly name
sg_object.assemblies[modified_input_bgc_name] = sg_object.assemblies.pop(input_bgc_name)
sg_object.assemblies[modified_input_bgc_name].id = modified_input_bgc_name

input_bgc_nuc_id = list(sg_object.assemblies[modified_input_bgc_name].loci.keys())[0]

# Annotate input BGC proteins with HMM models using external hmmscan
sg_object.annotate(
    hmm_directory=hmm_dir,
    use_neo4j_precalc=True,
)


sg_object.proteins.keys()

a = CompareDomains()

a.compare_all_to_all_parallel()
