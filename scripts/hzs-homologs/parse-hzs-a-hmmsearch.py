"""
Script to extract protein accessions with hit to HZS-A middle domain Pfam.
"""
import sys

# Parse output from hmmsearch
proteins = []
summary_data = {"assembly_accession":[], "rpoB_identifier":[], "rpoB_hmm_evalue":[]}
with open(sys.argv[1], "r") as f:
    for line in f:
        if line[0] != "#":
            line = line.strip().split()
            if float(line[4]) <= 1e-10:  # Remove hits with e-value less than 1e-30
                proteins.append(line[0])

print("hzs_a_middle_domain_containing_protein")
for protein in proteins:
    print(protein)
