"""
Script to extract taxa in the class Brocadiae in GTDB.
"""
import pandas as pd
import sys


# Read GTDB metadata
df = pd.read_csv(sys.argv[1], sep="\t")

# Split GTDB taxonomy and clean
df[
    ["domain", "phylum", "class", "order", "family", "genus", "species"]
] = df.gtdb_taxonomy.str.split(";", expand=True)

df["domain"] = df["domain"].str.replace("d__", "")
df["phylum"] = df["phylum"].str.replace("p__", "")
df["class"] = df["class"].str.replace("c__", "")
df["order"] = df["order"].str.replace("o__", "")
df["family"] = df["family"].str.replace("f__", "")
df["genus"] = df["genus"].str.replace("g__", "")
df["species"] = df["species"].str.replace("s__", "")

# Extract Brocadiae
df = df[df["class"] == 'Brocadiae']

# Write to file
df.to_csv(sys.argv[2], sep="\t", index=False)
