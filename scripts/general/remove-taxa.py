"""
Script to remove sequences from fasta file based on accession.
"""
from Bio import SeqIO
import sys

remove = []
with open(sys.argv[2]) as f:
    for line in f:
        remove.append(line.strip())
print(remove)
keep = []
for record in SeqIO.parse(sys.argv[1], "fasta"):
    if not record.id in remove:
        keep.append(record)

SeqIO.write(keep, sys.argv[3], "fasta")
