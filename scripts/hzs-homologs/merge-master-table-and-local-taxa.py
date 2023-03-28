import pandas as pd
import sys
from Bio import SeqIO

master_df = pd.read_csv(sys.argv[1], sep="\t")
local_taxa_df = pd.read_csv(sys.argv[2], sep="\t")

local_seq = []
for record in SeqIO.parse(sys.argv[3], "fasta"):
    local_seq.append(record.id)

local_taxa_df = local_taxa_df[local_taxa_df["taxa"].isin(local_seq)]
local_taxa_df.rename(columns={"stitle":"annotation"}, inplace=True)
master_df = pd.concat([master_df, local_taxa_df]).fillna("None")
master_df["tree_annotation"] = master_df["species"] + " - " + master_df["taxa"]
master_df.to_csv(sys.argv[4], sep="\t", index=False)

