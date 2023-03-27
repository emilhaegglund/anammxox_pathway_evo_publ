import pandas as pd
import seaborn as sns
from Bio import SeqIO
import sys
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

new_rc_params = {"text.usetex": False, "svg.fonttype": "none"}
mpl.rcParams.update(new_rc_params)

y_label_order = []

# Read file that contains the order of the assembly asseccions
# according to their position in the species phylogeny.
with open(sys.argv[1], "r") as f:
    for line in f:
        line = line.strip()
        y_label_order.append(line)
print(y_label_order)
query_length = {}

for record in SeqIO.parse(
    sys.argv[2],
    "fasta",
):
    query_length[record.id] = len(record.seq)
protein_df = pd.read_csv(sys.argv[3], sep="\t")
blast_header = [
    "query",
    "subject",
    "pid",
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
blast_df = pd.read_csv(
    sys.argv[4],
    names=blast_header,
    sep="\t",
)
blast_df.sort_values(by='bitscore', inplace=True, ascending=False)
blast_df.drop_duplicates(subset=["subject"], inplace=True)
anammox_pathway_proteins_df = pd.read_csv(
    sys.argv[5], sep="\t"
)
blast_df = pd.merge(
    left=blast_df, right=protein_df, left_on="subject", right_on="protein_id"
)
blast_df = pd.merge(
    left=blast_df,
    right=anammox_pathway_proteins_df,
    left_on="query",
    right_on="protein_id",
)
hit_map = {}
for accession, df in blast_df.groupby("accession"):
    accession = accession.strip()
    hit_map[accession] = {}
    for protein in anammox_pathway_proteins_df["protein_id"].to_list():
        hit_map[accession][protein] = 0
        if protein in df["query"].to_list():
            hits_df = df[df["query"] == protein].copy(deep=True)
            hits_df.sort_values(by="bitscore", inplace=True, ascending=False)
            best_hit = hits_df["bitscore"].iat[0]
            hit_map[accession][protein] = best_hit / query_length[protein]

for accession in y_label_order:
    if accession not in hit_map.keys():
        hit_map[accession] = {}
        for protein in anammox_pathway_proteins_df["protein_id"].to_list():
            hit_map[accession][protein] = 0

rename_index = {}
genome_df = pd.read_csv(sys.argv[6], sep="\t")
for accession in y_label_order:
    species_name = genome_df[genome_df["taxa"] == accession][
        "organism-name-modifier"
    ].iat[0]
    rename_index[accession] = species_name


df = pd.DataFrame.from_dict(hit_map)
df = df.T
df.replace(0, np.nan, inplace=True)
df = df.reindex(index=y_label_order)
rename_columns = {}
for protein_id in list(df.columns.values):

    func_annotation = anammox_pathway_proteins_df[
        anammox_pathway_proteins_df["protein_id"] == protein_id
    ]["functional_annotation"].iat[0]
    if not pd.isnull(func_annotation):
        rename_columns[protein_id] = func_annotation
    else:
        locus_tag = anammox_pathway_proteins_df[
            anammox_pathway_proteins_df["protein_id"] == protein_id
        ]["locus_tag"].iat[0]
        rename_columns[protein_id] = locus_tag
df.rename(index=rename_index, inplace=True)

plot_dfs = []
groups = [
    "HZS",
    "HAO-like",
    "Nir_associated",
    "NXR",
    #"ETM",
    #"Rb-2",
    #"Rb-3",
    "F-ATPase",
    #"AMT",
    #"FocA",
    #"NarK",
]
for group in groups:
    group_protein_id = anammox_pathway_proteins_df[
        anammox_pathway_proteins_df["group"] == group
    ]["protein_id"].to_list()
    plot_df = df[group_protein_id]
    plot_dfs.append(plot_df)

vmax = df.values.max()
vmin = df.values.min()

# Set the width for each group to the same size as the number of
# orthogroups in the group
plot_widths = [f.shape[1] for f in plot_dfs] + [1]
sns.set(rc={"figure.figsize": (5, 6)})
sns.set(font_scale=0.5)

fig, axs = plt.subplots(
    ncols=len(plot_dfs) + 1, gridspec_kw=dict(width_ratios=plot_widths)
)

cmap = sns.cm.rocket_r
for i, plot_df in enumerate(plot_dfs):
    plot_df.rename(columns=rename_columns, inplace=True)
    if i == 0:
        sns.heatmap(
            data=plot_df,
            square=True,
            cmap=cmap,
            ax=axs[i],
            cbar=False,
            linewidths=1,
        )
    else:
        sns.heatmap(
            data=plot_df,
            square=True,
            cmap=cmap,
            ax=axs[i],
            cbar=False,
            linewidths=1,
            yticklabels=False,
        )
    axs[i].tick_params(
        axis="both",
        which="major",
        labelsize=10,
        labelbottom=False,
        bottom=False,
        top=False,
        labeltop=True,
    )
    axs[i].set(xlabel=groups[i])
    axs[i].xaxis.set_label_position("top")
    # plot = sns.heatmap(data=df, cmap=cmap, square=True)
    # plt.tick_params(
    #    axis="both",
    #    which="major",
    #    labelsize=10,
    #    labelbottom=False,
    #    bottom=False,
    #    top=False,
    #    labeltop=True,
    # )

    xlabels = axs[i].get_xticklabels()
    axs[i].set_xticklabels(xlabels, rotation=-90)
fig.colorbar(axs[-2].collections[0], cax=axs[-1], shrink=0.5)
plt.savefig(sys.argv[7], bbox_inches="tight")
