import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys

hzs_b_df = pd.read_csv(sys.argv[1], sep="\t")
hzs_b_df["protein"] = "HZS-B"

hzs_c_df = pd.read_csv(sys.argv[2], sep="\t")
hzs_c_df["protein"] = "HZS-C"

hzs_a_df = pd.read_csv(sys.argv[3], sep="\t")
hzs_a_df["protein"] = "HZS-A"


hzs_df = pd.concat([hzs_b_df, hzs_c_df, hzs_a_df])
print(hzs_df.columns)
hzs_df = hzs_df[hzs_df["dS"] < 3]
hzs_df = hzs_df[hzs_df["dN"] > 0.01]
hzs_df = hzs_df[hzs_df["dnds"] < 10]

ax = sns.violinplot(data=hzs_df, x="protein", y="dnds")
ax.set(ylabel='dN/dS', xlabel="")

plt.savefig(sys.argv[4])
