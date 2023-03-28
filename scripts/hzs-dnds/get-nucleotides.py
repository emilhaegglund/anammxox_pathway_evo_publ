import sys
from Bio import SeqIO
import os

aa_headers = []
for record in SeqIO.parse(sys.argv[1], "fasta"):
    aa_headers.append(record.id)

store_record = []
for f in os.listdir(sys.argv[2]):
    f_path = os.path.join(sys.argv[2], f)
    for record in SeqIO.parse(f_path, "fasta"):
        if record.id.split("_")[2] in aa_headers:
            record.id = record.id.split("_")[2]
            store_record.append(record)

SeqIO.write(store_record, sys.argv[3], 'fasta')
