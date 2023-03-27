import pandas as pd
from Bio import SeqIO
import os
import sys

protein_accession_map = {"protein_id":[], "accession":[]}
for f in os.listdir(sys.argv[1]):
    accession = os.path.splitext(f)[0]
    f_path = os.path.join(sys.argv[1], f)
    for record in SeqIO.parse(f_path, "fasta"):
        protein_accession_map["protein_id"].append(record.id)
        protein_accession_map["accession"].append(accession)

for f in os.listdir(sys.argv[2]):
    accession = "_".join(f.split("_")[:2])
    f_path = os.path.join(sys.argv[2], f)
    for record in SeqIO.parse(f_path, "fasta"):
        protein_accession_map["protein_id"].append(record.id)
        protein_accession_map["accession"].append(accession)
df = pd.DataFrame.from_dict(protein_accession_map)

df.to_csv(sys.argv[3], sep="\t", index=False)

