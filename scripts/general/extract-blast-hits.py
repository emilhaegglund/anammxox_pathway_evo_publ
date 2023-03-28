import pandas as pd
import argparse
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--blast', type=str, required=True)
    parser.add_argument('--queries', type=str)
    parser.add_argument('--fasta', type=str, required=True)
    parser.add_argument('--anammox', action='store_true', default=False)
    parser.add_argument('--qcover', type=float),
    parser.add_argument('--scover', type=float),
    parser.add_argument('--random', type=int),
    parser.add_argument('--top', type=int),
    parser.add_argument('--source', type=str)
    return parser.parse_args()

args = parse_args()

if args.top and args.random:
    sys.exit("--top and --random are mutually exclusive")

blast_df = pd.read_csv(args.blast, sep='\t', compression="gzip")
print(blast_df.columns)
blast_df['qcover'] = blast_df["length"] / blast_df["qlen"]
blast_df['scover'] = blast_df["length"] / blast_df["slen"]

# Keep anammox
if not args.anammox:
    blast_df = blast_df[blast_df["class"] != "Brocadiae"]

# Select only entries from RefSeq
if args.source == "refseq":
    blast_df = blast_df[blast_df["sseqid"].str.contains("WP_")]

# Select only subset of queries
if args.queries:
    queries = args.queries.split(",")
    blast_df = blast_df[blast_df["qseqid"].isin(queries)]

# Filter sequences based on alignemnt length coverage
if args.qcover:
    blast_df = blast_df[blast_df["qcover"] >= args.qcover]
if args.scover:
    blast_df = blast_df[blast_df["scover"] >= args.scover]

# Remove duplicate entries
blast_df = blast_df.drop_duplicates(subset="sseqid", keep="first")

# Random sample N of blast hits
if args.random:
    random_dfs = []
    for query, query_df in blast_df.groupby("qseqid"):
        random_df = query_df.sample(args.random)
        random_dfs.append(random_df)
    blast_df = pd.concat(random_dfs)

# Select only N top blast hits
if args.top:
    top_dfs = []
    for query, query_df in blast_df.groupby("qseqid"):
        top_df = query_df.head(args.top)
        top_dfs.append(top_df)
    blast_df = pd.concat(top_dfs)

# Export fasta
with open(args.fasta, "w") as f:
    for i, row in blast_df.iterrows():
        f.write(f'>{row["stitle"]}\n')
        f.write(f'{row["full_sseq"].replace("*", "")}\n') 