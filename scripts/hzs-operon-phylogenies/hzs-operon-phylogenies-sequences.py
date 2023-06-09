import sys
import pandas as pd
from Bio import SeqIO


orthogroup = sys.argv[1]
orthogroup_df = pd.read_csv(sys.argv[2], sep="\t")
protein_ids = orthogroup_df[orthogroup_df["cluster_id"] == orthogroup][
    "feat_id"
].to_list()
if orthogroup == "OG0000174":
    protein_ids.append("GAX62881.1")
    protein_ids.remove("TWU54196.1")
store_records = []

for record in SeqIO.parse(sys.argv[3], 'fasta'):
    if record.id in protein_ids:
        store_records.append(record)

SeqIO.write(store_records, sys.argv[4], "fasta")

