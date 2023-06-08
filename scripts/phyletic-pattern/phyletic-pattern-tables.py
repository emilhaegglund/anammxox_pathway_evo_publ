import pandas as pd
import sys


blast_header = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
]
blast_df = pd.read_csv(sys.argv[1], sep="\t", names=blast_header)
protein_df = pd.read_csv(sys.argv[2], sep="\t")
genomes = []
with open(sys.argv[3], "r") as f:
    for line in f:
        genomes.append(line.strip())

blast_df = pd.merge(
    left=blast_df, right=protein_df, left_on="sseqid", right_on="protein_id"
)
print(blast_df[blast_df["sseqid"] == "SOH02641.1"])
query_seq = blast_df["qseqid"].unique().tolist()
data = {}
pident_data = {}
for query, query_df in blast_df.groupby("qseqid"):
    data[query] = {}
    pident_data[query + "_identity"] = {}
    query_df.sort_values(by="evalue", inplace=True)
    if query in ["SOH02641.1", "SOH02643.1"]:
        print(query_df[query_df["sseqid"] == "SOH02641.1"])
    for genome in genomes:
        if genome in query_df["accession"].to_list():
            hit = query_df[query_df["accession"] == genome]["sseqid"].iat[0]
            pident = query_df[query_df["accession"] == genome]["pident"].iat[0]
            data[query][genome] = hit
            pident_data[query + "_identity"][genome] = pident
        else:
            data[query][genome] = None
            pident_data[query + "_identity"][genome] = None

df = pd.DataFrame.from_dict(data)
pident_df = pd.DataFrame.from_dict(pident_data)
df = df.join(pident_df, how='inner')
print(df)
df_order = [
    "SOH05198.1",
    "SOH05199.1",
    "SOH05200.1",
    "SOH04660.1",
    "SOH04857.1",
    "SOH05157.1",
    "SOH03954.1",
    "SOH03957.1",
    "SOH03958.1",
    "SOH03959.1",
    "SOH02648.1",
    "SOH02643.1",
    "SOH02642.1",
    "SOH02641.1",
]
new_df_order = []
for i in df_order:
    new_df_order.append(i)
    new_df_order.append(i + "_identity")
df[new_df_order].to_csv(sys.argv[4], sep="\t")
