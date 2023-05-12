import pandas as pd
import sys

header = [
    "qseqid",
    "sseqid",
    "pident",
    "positive",
    "mismatch",
    "gapopen",
    "qlen",
    "slen",
    "length",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "stitle",
    "full_sseq",
]

df = pd.read_csv(sys.argv[1], sep="\t", names=header, compression="gzip")
queries = df["qseqid"].unique().tolist()
df = df[df["sseqid"].isin(queries)]

print(df[["qseqid", "sseqid", "pident"]])
identity_df = df.pivot(index="qseqid", columns="sseqid", values="pident")
identity_df.to_csv(sys.argv[2], sep="\t")
