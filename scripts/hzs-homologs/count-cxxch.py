import re
from Bio import SeqIO
import sys
import pandas as pd

p = re.compile("C..C[K,H]")
data = []
for record in SeqIO.parse(sys.argv[1], 'fasta'):
    counter = 0
    heme_pos = []
    for m in p.finditer(str(record.seq)):
        counter += 1
        heme_pos.append((m.start(), m.end()))
    data.append([record.id, heme_pos, counter])

heme_df = pd.DataFrame(data, columns=["protein_accession", "heme_pos", "n_heme_motifs"])

master_df = pd.read_csv(sys.argv[2], sep="\t")
master_df = pd.merge(master_df, heme_df, left_on="taxa", right_on="protein_accession", how="left")

master_df.to_csv(sys.argv[3], sep="\t", index=False)
