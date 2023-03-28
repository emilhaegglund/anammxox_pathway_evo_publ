import pandas as pd
import seaborn as sns
import sys
import matplotlib.pyplot as plt
from Bio import SeqIO
import os

sns.set(rc={"figure.figsize": (12, 12)})

# Import nucleotide identity data
df = pd.read_csv(sys.argv[1], sep="\t")


# The input data only contain each pair once, duplicate the dataframe and
# reverse the protein name columns and merge them to get all against all
df.rename(
    columns={"protein_id_i": "Protein_1", "protein_id_j": "Protein_2"}, inplace=True
)

# Create a matrix from the columns
df_wide = df.pivot_table(
    index="Protein_1", columns="Protein_2", values="pident", dropna=False
)
proteins = df_wide.index.to_list()

# Read genome info
genome_df = pd.read_csv(sys.argv[3], sep="\t")

# Sort based on the order in th rpoB-phylogeny
sort_order = []
y_labels = []
sort_accession_order = [
    "GCA_002443295.1_ASM244329v1_protein.faa",
    "GCA_013112645.1_ASM1311264v1_protein.faa",
    "GCA_900232105.1_Kuenenia_stuttgartiensis_MBR1_protein.faa",
    "GCA_011066545.1_ASM1106654v1_protein.faa",
    "GCA_021650915.1_ASM2165091v1_protein.faa",
    "GCA_017347445.1_ASM1734744v1_protein.faa",
    "GCA_000296795.1_ASM29679v1_protein.faa",
    "GCA_021650895.1_ASM2165089v1_protein.faa",
    "GCA_000949635.1_ASM94963v1_protein.faa",
    "GCA_007860005.1_ASM786000v1_protein.faa",
]
for f in sort_accession_order:
    accession = "_".join(f.split("_")[:2])
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

# Change the labels on the y-axis
g.set_yticklabels(y_labels)
plt.savefig(sys.argv[4], bbox_inches="tight")
