"""
Script to add column with locus tags for the blast results
between MBR1 and kust.
"""
from Bio import SeqIO
import pandas as pd
import sys


# Create a dataframe with mapping between protein id and locus tags
data = {"protein_id":[], "kust_locus_tag":[]}
for record in SeqIO.parse(sys.argv[2], "genbank"):
    final_features = []
    for item in record.features:
        if item.type == 'CDS':
            locus_tag = item.qualifiers['locus_tag'][0]
            protein_id = item.qualifiers['protein_id'][0]
            data["protein_id"].append(protein_id)
            data["kust_locus_tag"].append(locus_tag)

locus_tag_df = pd.DataFrame.from_dict(data)

# Read blast
blast_header = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
]
blast_df = pd.read_csv(sys.argv[1], sep="\t", names=blast_header)

# Merge the dataframes and save
blast_df = pd.merge(blast_df, locus_tag_df, left_on="sseqid", right_on="protein_id")
blast_df.drop(columns="protein_id", inplace=True)
blast_df.to_csv(sys.argv[3], sep="\t", index=False)
