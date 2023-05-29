import pandas as pd
import sys

blast_headers = [
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
    "full_sseq"
]
df = pd.read_csv(sys.argv[1], sep='\t', compression="gzip", names=blast_headers)
print(df.shape)
df.drop_duplicates(subset=["sseqid"], inplace=True)

for i, row in df.iterrows():
    print(">" + row["stitle"])
    print(row["full_sseq"].replace("*", ""))
