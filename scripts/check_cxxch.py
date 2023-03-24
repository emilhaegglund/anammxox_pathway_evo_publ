import re
from Bio import SeqIO
import sys
import pandas as pd

p = re.compile("C..C[H,K]")
data = []
store_records = []
for record in SeqIO.parse(sys.argv[1], 'fasta'):
    counter = 0
    heme_pos = []
    for m in p.finditer(str(record.seq)):
        counter += 1
        heme_pos.append((m.start(), m.end()))
    data.append([record.id, heme_pos, counter])
    if counter == 8:
        store_records.append(record)

heme_df = pd.DataFrame(data, columns=["protein_accession", "heme_pos", "n_heme_motifs"])

heme_df.to_csv(sys.argv[2], sep="\t", index=False)
SeqIO.write(store_records, sys.argv[3], "fasta")
