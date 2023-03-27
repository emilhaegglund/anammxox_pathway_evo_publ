"""
Script to extract the protein with the best hit to rpoB for each genome.
"""
from Bio import SeqIO
import sys
import os
import pandas as pd


# Read protein mappings file
df = pd.read_csv(sys.argv[1], sep="\t")

# Read all proteins
protein_dict = SeqIO.to_dict(SeqIO.parse(sys.argv[2], "fasta"))

# Parse output from hmmsearch
species_dict = {}
proteins = []
summary_data = {"assembly_accession":[], "rpoB_identifier":[], "rpoB_hmm_evalue":[]}
with open(sys.argv[3], "r") as f:
    for line in f:
        if line[0] != "#":
            line = line.strip().split()
            if float(line[4]) <= 1e-30:  # Remove hits with e-value less than 1e-30
                assembly_accession = df[df["protein_accession"] == line[0]][
                    "assembly_accession"
                ].iat[0]
                # Store the best hit
                if assembly_accession not in species_dict.keys():
                    species_dict[assembly_accession] = line[0]
                    proteins.append(line[0])
                    summary_data['assembly_accession'].append(assembly_accession)
                    summary_data['rpoB_identifier'].append(line[0])
                    summary_data['rpoB_hmm_evalue'].append(float(line[4]))
                    

# Extract proteins and write to file
rpoB_records = []
for assembly_accession, protein_accession in species_dict.items():
    record = protein_dict[protein_accession]
    record.id = assembly_accession
    rpoB_records.append(record)
SeqIO.write(rpoB_records, sys.argv[4], "fasta")

summary_df = pd.DataFrame.from_dict(summary_data)
summary_df.to_csv(sys.argv[5], sep="\t", index=False)
