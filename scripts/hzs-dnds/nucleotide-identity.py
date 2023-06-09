from Bio import SeqIO
import sys
import pandas as pd


protein_ids = []
alignments = []
for record in SeqIO.parse(sys.argv[1], "fasta"):
    protein_ids.append(record.id)
    alignments.append(record.seq)

data = {"protein_id_i":[], "protein_id_j":[], "identical":[], "pident":[]}
for i, protein_id_i in enumerate(protein_ids):
    for j, protein_id_j in enumerate(protein_ids):
        count = 0
        aln_length = 0
        for pos in range(len(alignments[i])):
            nucl_i = alignments[i][pos]
            nucl_j = alignments[j][pos]
            if nucl_i != "-" and nucl_j != "-":
                if nucl_i == nucl_j:
                    count += 1
                aln_length += 1
        data["protein_id_i"].append(protein_id_i)
        data["protein_id_j"].append(protein_id_j)
        data["identical"].append(count)
        if aln_length > 0:
            data["pident"].append(count/aln_length)
        else:
            data["pident"].append(0)


df = pd.DataFrame.from_dict(data)
df.to_csv(sys.argv[2], sep="\t", index=False)

