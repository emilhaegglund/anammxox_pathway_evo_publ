import sys


seq_to_extract = sys.argv[1]
clusters_path = sys.argv[2]

repr_seqs = []
with open(seq_to_extract, "r") as f:
    for line in f:
        line = line.strip()
        repr_seqs.append(line)

cluster = []
seqs_to_keep = []
with open(clusters_path, "r") as f:
    for line in f:
        line = line.strip()
        if line[0] == ">":
            for seq in repr_seqs:
                if seq in cluster:
                    seqs_to_keep += cluster
            cluster = []
        else:
            seq = line.split()[2]
            seq = seq.replace(">", "")
            seq = seq.replace("...", "")
            cluster.append(seq)

for seq in seqs_to_keep:
    print(seq)
