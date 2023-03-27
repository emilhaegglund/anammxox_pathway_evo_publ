import pandas as pd
import sys
import os



assemblies_df = pd.read_csv(sys.argv[1], sep="\t")
print(assemblies_df.columns)
col_rename = {"Assembly Accession":"assembly_accession",
              "Organism Name":"ncbi_organism_name",
              "Organism Infraspecific Names Strain":"ncbi_strain",
              "Organism Infraspecific Names Isolate":"ncbi_isolate",
              "Assembly Stats GC Percent":"ncbi_gc_content",
              "Assembly Level":"ncbi_assembly_level",
              "Assembly Stats Number of Contigs":"ncbi_contig_count",
              "Assembly Stats Total Sequence Length":"ncbi_total_sequence_length"}
assemblies_df.rename(columns=col_rename, inplace=True)

# Read NCBI
ncbi_df = pd.read_csv(sys.argv[2], sep="\t")
ncbi_taxa = ncbi_df["Assembly Accession"].to_list()
assemblies_df["ncbi_taxonomy"] = assemblies_df["assembly_accession"].map(
    lambda x: x in ncbi_taxa
)

# Read GTDB
gtdb_df = pd.read_csv(sys.argv[3], sep="\t")
gtdb_taxa = gtdb_df["ncbi_genbank_assembly_accession"].to_list()
assemblies_df["gtdb_taxonomy"] = assemblies_df["assembly_accession"].map(
    lambda x: x in gtdb_taxa
)

# Read Prokka
prokka_annotation = []
for p in os.listdir(sys.argv[4]):
    prokka_annotation.append(p)
assemblies_df["prokka_annotation"] = assemblies_df["assembly_accession"].map(
    lambda x: x in prokka_annotation
)

rpoB_df = pd.read_csv(sys.argv[5], sep="\t")
assemblies_df = pd.merge(
    left=assemblies_df, right=rpoB_df, on="assembly_accession", how="left"
)

pvc = []
with open(sys.argv[6], "r") as f:
    for line in f:
        pvc.append(line.strip())

extended_dataset = []
with open(sys.argv[7], "r") as f:
    for line in f:
        extended_dataset.append(line.strip())

extended_dataset += pvc
assemblies_df["dataset_extended"] = assemblies_df["assembly_accession"].map(
    lambda x: x in extended_dataset
) 

hq_dataset = []
with open(sys.argv[8], "r") as f:
    for line in f:
        hq_dataset.append(line.strip())

with open(sys.argv[9], "r") as f:
    for line in f:
        hq_dataset.append(line.strip())
assemblies_df["dataset_high_quality"] = assemblies_df["assembly_accession"].map(
    lambda x: x in hq_dataset
) 
selected_cols = ["ncbi_genbank_assembly_accession", "checkm_completeness", "checkm_contamination"]
assemblies_df["outgroup"] = assemblies_df["assembly_accession"].map(lambda x: x in pvc)
assemblies_df = pd.merge(left=assemblies_df, right=gtdb_df[selected_cols], left_on="assembly_accession", right_on="ncbi_genbank_assembly_accession", how="left")
assemblies_df.rename(columns={"checkm_completeness":"gtdb_checkm_completeness", "checkm_contamination":"gtdb_checkm_contamination"}, inplace=True)
assemblies_df.drop(columns=["ncbi_genbank_assembly_accession"], inplace=True)
assemblies_df.to_csv(
    "/media/argos-emiha442/emiha442/tmp/rpoB_table_1.tsv", index=False, sep="\t"
)
