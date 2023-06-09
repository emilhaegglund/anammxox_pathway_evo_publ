import pandas as pd
import seaborn as sns
import sys
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
import os

sns.set(rc={"figure.figsize": (12, 12)})

# Import dnds data
df = pd.read_csv(sys.argv[1], sep="\t")

# Exclude all dnds-values that has a dS larger than 3 and dN less than 0.01
df.loc[(df["dS"] > 3) | (df["dN"] < 0.01), "dnds"] = np.nan

# The input data only contain each pair once, duplicate the dataframe and
# reverse the protein name columns and merge them to get all against all
df_2 = df.copy()
df_2.rename(columns={"Gene_1": "Gene_2", "Gene_2": "Gene_1"}, inplace=True)
df = pd.concat([df, df_2])

# Create a matrix from the columns
df_wide = df.pivot_table(index="Gene_1", columns="Gene_2", values="dnds", dropna=False)
proteins = df_wide.index.to_list()
# Sort based on the order in th rpoB-phylogeny
genome_df = pd.read_csv(sys.argv[3], sep="\t")
sort_order = []
y_labels = []
sort_accession_order = [
    "GCA_002443295.1_ASM244329v1_protein.faa",
    "GCA_013112645.1_ASM1311264v1_protein.faa",
    "GCA_900232105.1_Kuenenia_stuttgartiensis_MBR1_protein.faa",
    "GCA_021650915.1_ASM2165091v1_protein.faa",
    "GCA_017347445.1_ASM1734744v1_protein.faa",
    "GCA_000296795.1_ASM29679v1_protein.faa",
    "GCA_021650895.1_ASM2165089v1_protein.faa",
    "GCA_000949635.1_ASM94963v1_protein.faa",
    "GCA_007860005.1_ASM786000v1_protein.faa",
]
for f in sort_accession_order:
    accession = "_".join(f.split("_")[:2])  # Extract the genome accession from the file name
    f_path = os.path.join(sys.argv[2], f)
    for record in SeqIO.parse(f_path, "fasta"):
        if record.id in proteins:
            sort_order.append(record.id)
            publ_name = genome_df[genome_df["accession"] == accession]["publ_name"].iat[
                0
            ]
            y_labels.append(publ_name + " " + record.id)

df_wide = df_wide.reindex(sort_order)
df_wide = df_wide.reindex(sort_order, axis="columns")
g = sns.heatmap(
    df_wide, cmap="viridis", annot=True, linewidths=0.1, annot_kws={"fontsize": 8}
)
g.set_yticklabels(y_labels)
plt.savefig(sys.argv[4], bbox_inches="tight")
