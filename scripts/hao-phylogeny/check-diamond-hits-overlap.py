import pandas as pd
import sys

df = pd.read_csv(sys.argv[1], sep="\t", compression="gzip")
df = df[df["sseqid"].str.contains("WP_")]

df['qcover'] = df["length"] / df["qlen"]
df['scover'] = df["length"] / df["slen"]

df = df[df['qcover'] >= 0.5]
df = df[df['scover'] >= 0.5]

sim_dict = {}
for prot_1 in df["qseqid"].unique().tolist():
    for prot_2 in df["qseqid"].unique().tolist():
        prot_1_hits = set(df[df["qseqid"] == prot_1]["sseqid"].tolist())
        prot_2_hits = set(df[df["qseqid"] == prot_2]["sseqid"].tolist())
        ji = len(prot_1_hits.intersection(prot_2_hits))/len(prot_1_hits.union(prot_2_hits))
        print(prot_1, prot_2, len(prot_2_hits.intersection(prot_1_hits)))
        if prot_1 in sim_dict.keys():
            sim_dict[prot_1][prot_2] = ji
        else:
            sim_dict[prot_1] = {prot_2:ji}


sim_df = pd.DataFrame.from_dict(sim_dict)
sim_df.to_csv(sys.argv[2], sep="\t")


