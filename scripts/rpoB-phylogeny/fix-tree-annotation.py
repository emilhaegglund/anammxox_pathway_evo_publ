"""
Script to fix tree annotation file because of inconsistency in NCBI
"""
import pandas as pd
import sys

df = pd.read_csv(sys.argv[1], sep="\t")

data = {"taxa": [], "organism-name-modifier": []}
for i, row in df.iterrows():
    data["taxa"].append(row["Assembly Accession"])
    ending = row["Organism Name"].split()[-1]
    # If the last word of the organism name is same as strain or isolate,
    # keep as is. Otherwise merge organism name with strain if present.
    # If it's not present merge with isolate. Else just keep organism name as is.
    if ending == row["Organism Infraspecific Names Strain"]:
        data["organism-name-modifier"].append(row["Organism Name"])
    elif ending == row["Organism Infraspecific Names Isolate"]:
        data["organism-name-modifier"].append(row["Organism Name"])
    else:
        if type(row["Organism Infraspecific Names Strain"]) != float:
            data["organism-name-modifier"].append(
                row["Organism Name"] + " " + row["Organism Infraspecific Names Strain"]
            )
        elif type(row["Organism Infraspecific Names Isolate"]) != float:
            data["organism-name-modifier"].append(
                row["Organism Name"]
                + " "
                + row["Organism Infraspecific Names Isolate"].replace(".fa.gz", "")
            )
        else:
            data["organism-name-modifier"].append(row["Organism Name"])


annotation_df = pd.DataFrame.from_dict(data)
annotation_df.to_csv(sys.argv[2], sep="\t", index=False)
