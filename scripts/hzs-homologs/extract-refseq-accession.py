import pandas as pd
import sys

df = pd.read_csv(sys.argv[1], sep='\t')
if sys.argv[2] == "hzs_bc":
    df = df[df["qseqid"].isin(["GAX62881.1", "QII12198.1", "QII12199.1"])]
elif sys.argv[2] == "hzs_a":
    df = df[df["qseqid"].isin(["GAX62882.1", "QII12200.1"])]

df = df[df["ncbi_genome_category"] == "none"]

df = df.drop_duplicates(subset="assembly_accession", keep="first")

for protein in df["taxa"].tolist():
    accession = df[df["taxa"] == protein]["assembly_accession"].iloc[0]
    if protein[:3] != "WP_":
        accession = accession.replace("GCF","GCA")
    print(accession)