import pandas as pd
import sys

df = pd.read_csv(sys.argv[1], sep="\t")
df = df.drop_duplicates(subset="sseqid")

df.rename(columns={"sseqid":"taxa", "stitle":"annotation", "pident":"blast_pident", "evalue":"blast_evalue", "bitscore":"blast_bitscore"}, inplace=True)

output_columns = [
    "taxa",
    "qseqid",
    "annotation",
    "qlen",
    "slen",
    "blast_pident",
    "blast_evalue",
    "blast_bitscore",
    "domain",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
    "assembly_accession"
]

df[output_columns].to_csv(sys.argv[2], index=False, sep="\t")