import pandas as pd
import sys

df = pd.read_csv(sys.argv[1], sep='\t')
print(df.columns)
gtdb_df = pd.read_csv(sys.argv[2], sep="\t")

df["accession_wo_version"] = df["assembly_accession"].str.replace("GCF", "GCA")
df["accession_wo_version"] = df["accession_wo_version"].str.replace(r".[0-9]$", "")

gtdb_df["accession_wo_version"] = gtdb_df["ncbi_genbank_assembly_accession"].str.replace(r".[0-9]$", "")
gtdb_columns = ["accession_wo_version", "ncbi_genome_category"]
df = pd.merge(df, gtdb_df[gtdb_columns], on="accession_wo_version", how="left")


df.to_csv(sys.argv[3], sep="\t", index=False)