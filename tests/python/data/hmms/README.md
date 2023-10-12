conda install pfam_scan

conda install hmmer

wget <http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz>
wget <http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/active_site.dat.gz>
wget <http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz>

gunzip *

hmmpress Pfam-A.hmm

pfam_scan.pl -fasta <fasta_file> -dir <directory location of Pfam files>

to get: lagriamide_mibig_bgc0001946_pfam_online_annotations_2022-01-11.csv

<http://pfam.xfam.org> has some prreacaluclated PFAM annotations, including most?all? of UNIPROT

So I downloaded the associated JSON on <http://pfam.xfam.org> and parsed it with python:

import json
from pathlib import Path
import pandas as pd

temp_list = []

for filename in Path('/home/chase/Downloads/pfam_download').glob('*.json'):

# Opening JSON file

```python
with open(filename) as f:
    data = json.load(f)
    for i in data['regions']:
        temp_dict = i['metadata']
        temp_dict["protein"]= data['metadata']['accession']
        temp_list.append(temp_dict)

df = pd.DataFrame(temp_list)
df.to_csv("lagriamide_mibig_bgc0001946_pfam_online_annotations_2022-01-11.csv")
```

The HMM model's hash was then found and exported from the socialgene HMM file useing HMMER's `hmmfetch` command
