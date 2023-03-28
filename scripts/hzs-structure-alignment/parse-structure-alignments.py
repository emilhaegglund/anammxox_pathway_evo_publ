import sys
import os
import pandas as pd

data = {"reference":[], "target":[], "tm_score_ref":[], "tm_score_target":[], "tm_score_avg":[]}
for f_path in sys.argv[1:-2]:
    path = f_path.split(os.sep)
    target_structure = path[-1]
    ref_structure = path[-2].replace('_structure_alignment', '')
    print(ref_structure)
    data["reference"].append(ref_structure)
    data["target"].append(target_structure)
    tm_scores = []
    with open(f_path, "r") as f:
        for line in f:
            line = line.split()
            if len(line) > 1:
                if line[0] == "TM-score=":
                    tm_scores.append(float(line[1]))
    data["tm_score_ref"].append(tm_scores[0])
    data["tm_score_target"].append(tm_scores[1])
    data["tm_score_avg"].append(tm_scores[2])

df = pd.DataFrame.from_dict(data)
df.to_csv(sys.argv[-1], sep="\t", index=False)
