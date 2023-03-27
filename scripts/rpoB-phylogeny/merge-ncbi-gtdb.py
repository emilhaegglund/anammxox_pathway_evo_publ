"""
Script to merge accessions under Candidatus Brocadiia
in NCBI taxonomy and accessions under Brocadiae in GTDB taxonomy.
"""
import pandas as pd
import sys


# Read NCBI data
ncbi_df = pd.read_csv(sys.argv[1], sep="\t")
ncbi_assembly_accessions = ncbi_df["Assembly Accession"].to_list()
ncbi_assembly_accessions = set(ncbi_assembly_accessions)
print(f"ncbi_taxonomy anammox: {len(ncbi_assembly_accessions)}")

# Read GTDB data
gtdb_df = pd.read_csv(sys.argv[2], sep="\t")
gtdb_assembly_accessions = gtdb_df["ncbi_genbank_assembly_accession"].to_list()
gtdb_assembly_accessions = set(gtdb_assembly_accessions)
print(f"gtdb_taxonomy anammox: {len(gtdb_assembly_accessions)}")

# Combine taxa from NCBI taxonomy and GTDB taxonomy
all_anammox_accessions = ncbi_assembly_accessions.union(gtdb_assembly_accessions)
print(f"all anammox: {len(all_anammox_accessions)}")

# Write output to file
with open(sys.argv[3], "w") as f_out:
    for accession in all_anammox_accessions:
        f_out.write(accession + "\n")


