"""
Script to remove sequences from fasta file based on accession.
"""
from Bio import SeqIO
import sys

remove = ["GCA_905339135.1", "GCA_014384865.1"]  # sequences to remove
keep = []
for record in SeqIO.parse(sys.argv[1], "fasta"):
    if not record.id in remove:
        keep.append(record)

SeqIO.write(keep, sys.argv[2], "fasta")
