from Bio import SeqIO
import sys


seq_to_keep = []
for record in SeqIO.parse(sys.argv[1], "fasta"):
    if record.id != "GAX62881.1":
        seq_to_keep.append(record)

SeqIO.write(seq_to_keep, sys.argv[2], "fasta")
