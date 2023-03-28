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
blast_df = pd.read_csv(
    sys.argv[1],
    sep="\t",
    names=blast_headers,
    compression="gzip",
)
protein_df = pd.read_csv(sys.argv[2], sep="\t", compression="gzip")

blast_df = pd.merge(blast_df, protein_df, left_on="sseqid", right_on="protein_id")
blast_df.to_csv(
    sys.argv[3],
    sep="\t",
    compression="gzip",
    index=False,
)
