"""
Script to create the subset used in the analysis
"""
from Bio import SeqIO
import sys

# Read file with anammox accessions to keep for the subset
anammox_accessions = []
with open(sys.argv[1], "r") as f:
    for line in f:
        anammox_accessions.append(line.strip())

# Read file with PVC accessions to keep for the subset
pvc_accessions = []
with open(sys.argv[2], "r") as f:
    for line in f:
        pvc_accessions.append(line.strip())


outgroup_accessions = []
with open(sys.argv[3], "r") as f:
    for line in f:
        outgroup_accessions.append(line.strip())

keep_accessions = anammox_accessions + pvc_accessions + outgroup_accessions

# Filter original sequence file with RpoB sequences
keep = []
for record in SeqIO.parse(sys.argv[4], "fasta"):
    if record.id in keep_accessions:
        keep.append(record)

# Write subset of sequences to file
SeqIO.write(keep, sys.argv[5], "fasta")
