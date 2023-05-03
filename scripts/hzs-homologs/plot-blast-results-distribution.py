import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import sys

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
font = {'family' : 'Arial',
        'size'   :6}

matplotlib.rc('font', **font)


hzs_df = pd.read_csv(sys.argv[1], sep='\t', compression="gzip")
hzs_df = hzs_df[hzs_df["bitscore"] >= 80]
hzs_a_df = hzs_df[hzs_df["qseqid"].isin(["GAX62882.1", "SOH05200.1"])]
#hzs_a_df = hzs_a_df[hzs_a_df["sseqid"].str.contains("WP_")]
hzs_bc_df = hzs_df[hzs_df["qseqid"].isin(["GAX62881.1", "SOH05198.1", "SOH05199.1"])]
#hzs_bc_df = hzs_bc_df[hzs_bc_df["sseqid"].str.contains("WP_")]
hzs_a_df = hzs_a_df[hzs_a_df["class"] != "Brocadiae"]
hzs_bc_df = hzs_bc_df[hzs_bc_df["class"] != "Brocadiae"]
hzs_bc_fused_df = hzs_bc_df[hzs_bc_df["qseqid"] == "GAX62881.1"]
hzs_b_df = hzs_bc_df[hzs_bc_df["qseqid"] == "SOH05198.1"]
hzs_c_df = hzs_bc_df[hzs_bc_df["qseqid"] == "SOH05199.1"]
print(hzs_a_df.columns)

fig, axes = plt.subplots(nrows=3, ncols=4, figsize=(10, 6))
fig.suptitle('Summary BLAST results', fontsize=16)
print(axes)
sns.histplot(data=hzs_a_df, ax=axes[0,0], x="pident")
sns.histplot(data=hzs_b_df, ax=axes[0,1], x="pident")
sns.histplot(data=hzs_c_df, ax=axes[0,2], x="pident")
sns.histplot(data=hzs_bc_df, ax=axes[0,3], x="pident")
sns.histplot(data=hzs_a_df, ax=axes[1,0], x="bitscore")
sns.histplot(data=hzs_b_df, ax=axes[1,1], x="bitscore")
sns.histplot(data=hzs_c_df, ax=axes[1,2], x="bitscore")
sns.histplot(data=hzs_bc_df, ax=axes[1,3], x="bitscore")
sns.histplot(data=hzs_a_df, ax=axes[2,0], x="length")
sns.histplot(data=hzs_b_df, ax=axes[2,1], x="length")
sns.histplot(data=hzs_c_df, ax=axes[2,2], x="length")
sns.histplot(data=hzs_bc_df, ax=axes[2,3], x="length")
axes[0,0].set_title("HZS-A (807AA)")
axes[0,1].set_title("HZS-B (377AA)")
axes[0,2].set_title("HZS-C (341AA)")
axes[0,3].set_title("HZS-BC (690AA)")
for i in range(0,3):
    for j in range(1,4):
        print(i,j)
        axes[i,j].set_ylabel("")
plt.tight_layout()
plt.savefig(sys.argv[2])
