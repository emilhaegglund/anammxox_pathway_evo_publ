from Bio import SeqIO, SeqRecord
import sys


seq_to_keep = []
for record in SeqIO.parse(sys.argv[1], "fasta"):
    if record.id != "TWU54196.1":
        seq_to_keep.append(record)

SeqIO.write(seq_to_keep, sys.argv[2], "fasta")
