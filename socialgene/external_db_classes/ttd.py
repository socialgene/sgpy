from socialgene.base.socialgene import SocialGene


import re

sg_obj = SocialGene()
fapath = "/home/chase/Downloads/ttd_database/P2-06-TTD_sequence_all.txt"

with open(fapath, "rt") as h:
    sequence = ""
    after_header = False
    for line in h:
        if line.startswith(">"):
            after_header = True
            if sequence:
                sg_obj.add_protein(
                    sequence=sequence,
                    external_protein_id=protein_id,
                )
            protein_id = re.search("T[0-9]{5}", line).group()
            sequence = ""
        elif after_header:
            sequence += line.strip()
