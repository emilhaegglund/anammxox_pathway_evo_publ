import pandas as pd
import sys
import seaborn as sns
import matplotlib.pyplot as plt


# Read data and prepare for plot
df = pd.read_csv(sys.argv[1], sep="\t")

# Create a best column with the best score from the ref, target and avg columns
df["tm_score_best"] = df[["tm_score_ref", "tm_score_target", "tm_score_avg"]].max(axis=1) 
df_plot = pd.melt(df, id_vars=["reference", "target"],
       value_vars=["tm_score_ref",
                   "tm_score_target",
                   "tm_score_avg",
                   "tm_score_best",
                   ],
       var_name='tm_method',
       value_name="tm_score")

# Plot boxplots of the TM-scores for reference protein length,
# target protein length, and average protein length.
sns_plot = sns.boxplot(data=df_plot, x="tm_method", y="tm_score", hue="reference")
sns_plot.set(xlabel='Method',
       ylabel='TM-Score')
plt.savefig(sys.argv[2])


