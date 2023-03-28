import pandas as pd
import sys

blast_path = sys.argv[1]
gtdb_path = sys.argv[2]
output_path = sys.argv[3]

# Read blast table
blast_df = pd.read_csv(blast_path, sep="\t", compression="gzip")
print(blast_df.shape)

# Read GTDB-metadata
df = pd.read_csv(gtdb_path, sep="\t")

# Subset the metadata
df = df[
    [
        "assembly_accession",
        "domain",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
    ]
]

# Fixing RefSeq/Genbank assembly accessions due to the cleaning before downloading GTDB-data
blast_accessions = blast_df["assembly_accession"].to_list()
all_assembly_accessions = df["assembly_accession"].to_list()
missing_accessions = set(blast_accessions) - set(all_assembly_accessions)
m = blast_df["assembly_accession"].isin(list(missing_accessions))
blast_df.loc[m, "assembly_accession"] = blast_df.loc[
    m, "assembly_accession"
].str.replace("GCA", "GCF")

# Add taxonomy to blast results
df = pd.merge(blast_df, df, on="assembly_accession")
print(df.shape)
df.to_csv(output_path, sep="\t", compression="gzip", index=False)
