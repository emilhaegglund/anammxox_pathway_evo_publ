"""
Script to find AlphaFold Structurs for homologs to HZS.
"""
import requests
import json
import pandas as pd
import sys

# Helper function to download data
# https://colab.research.google.com/drive/1i9UtVqa4m9WQ4ZVJkbGdwWbto_W7zmP_#scrollTo=xMmq9GVxzmWl
def get_url(url, **kwargs):
    response = requests.get(url, **kwargs)

    if not response.ok:
        print(response.text)
        response.raise_for_status()
        sys.exit()

    return response


WEBSITE_API = "https://rest.uniprot.org/"

df = pd.read_csv(
    sys.argv[1],
    sep="\t",
    compression="gzip",
)
if sys.argv[3] == "hzs_bc":
    df = df[df["qseqid"].isin(["GAX62881.1", "SOH05198.1", "SOH05199.1"])]
elif sys.argv[3] == "hzs_a":
    df = df[df["qseqid"].isin(["GAX62882.1", "SOH05200.1"])]
else:
    sys.exit()

# Remove bad hits
df = df[df["bitscore"] >= 80]

# Dictionary to store results from Uniprot
results = {"accession": [], "alphafold_url": []}

# For each hit in the blast-search, see if
# the protein sequence have a AlphaFold-predicted structure
for hit in df["sseqid"]:
    r = get_url(f"{WEBSITE_API}/uniprotkb/search?query={hit}")
    data = r.json()
    results["accession"].append(hit)
    if len(data["results"]) > 0:
        inAlphaFold = False
        for db in data["results"][0]["uniProtKBCrossReferences"]:
            if db["database"] == "AlphaFoldDB":
                accession = data["results"][0]["primaryAccession"]
                af_url = (
                    f"https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v4.pdb"
                )
                print(
                    f"{hit}\thttps://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v4.pdb"
                )
                results["alphafold_url"].append(af_url)
                inAlphaFold = True
        if not inAlphaFold:  # If there is no predicted structure
            print(f"{hit} not in AlphaFoldDB")
            results["alphafold_url"].append(None)
    else:  # If hit is not in UniProt
        print(f"{hit} not in Uniprot")
        results["alphafold_url"].append(None)

# Write results to file
results_df = pd.DataFrame.from_dict(results)
results_df.to_csv(sys.argv[2], sep="\t", index=False)
