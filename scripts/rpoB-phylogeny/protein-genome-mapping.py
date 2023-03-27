"""
Script to create a file that maps protein accession to assembly accession.
This is needed to map the RpoB proteins to the correct taxa.
"""
import sys
import os
import pandas as pd
from Bio import SeqIO

# Loop over all proteomes in the input
protein_mapping = {"protein_accession": [], "assembly_accession": []}
for f in sys.argv[1:-1]:
    proteome = os.path.split(f)[1]
    # Check if the proteome is from NCBI or from Prokka
    if proteome == "protein.faa":  # NCBI proteomes
        accession = os.path.dirname(f)
        accession = os.path.split(accession)[1]
        print(accession)
    else:  # Prokka proteomes
        accession = "_".join(proteome.split("_")[:2])
        print(accession)
    for record in SeqIO.parse(f, "fasta"):
        protein_mapping["protein_accession"].append(record.id)
        protein_mapping["assembly_accession"].append(accession)

# Convert to dataframe and save as tsv
df = pd.DataFrame.from_dict(protein_mapping)
df.to_csv(sys.argv[-1], sep="\t", index=False)
