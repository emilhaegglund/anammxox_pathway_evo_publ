import pandas as pd
from Bio import SeqIO
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
              "Assembly Stats Total Sequence Length":"ncbi_total_genome_size"}
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

# Read RpoB summary
rpoB_df = pd.read_csv(sys.argv[5], sep="\t")
assemblies_df = pd.merge(
    left=assemblies_df, right=rpoB_df, on="assembly_accession", how="left"
)

# Read PVC-accessions
pvc = []
with open(sys.argv[6], "r") as f:
    for line in f:
        pvc.append(line.strip())

# Read extended dataset
extended_dataset = []
with open(sys.argv[7], "r") as f:
    for line in f:
        extended_dataset.append(line.strip())

extended_dataset += pvc
assemblies_df["rpoB_extended"] = assemblies_df["assembly_accession"].map(
    lambda x: x in extended_dataset
)

# Read HQ-dataset
hq_dataset = []
with open(sys.argv[8], "r") as f:
    for line in f:
        hq_dataset.append(line.strip())

# Read PVC subset
with open(sys.argv[9], "r") as f:
    for line in f:
        hq_dataset.append(line.strip())
assemblies_df["rpoB_high_quality"] = assemblies_df["assembly_accession"].map(
    lambda x: x in hq_dataset
)

rpoB_all = []
for record in SeqIO.parse(sys.argv[10], "fasta"):
    rpoB_all.append(record.id)
assemblies_df["rpoB_all"] = assemblies_df["assembly_accession"].isin(rpoB_all)

rpoB_cleaned = []
for record in SeqIO.parse(sys.argv[11], "fasta"):
    rpoB_cleaned.append(record.id)
assemblies_df["rpoB_wo_sporious_taxa"] = assemblies_df["assembly_accession"].isin(rpoB_cleaned)

selected_cols = ["ncbi_genbank_assembly_accession", "checkm_completeness", "checkm_contamination"]
assemblies_df["outgroup"] = assemblies_df["assembly_accession"].map(lambda x: x in pvc)
assemblies_df = pd.merge(left=assemblies_df, right=gtdb_df[selected_cols], left_on="assembly_accession", right_on="ncbi_genbank_assembly_accession", how="left")
assemblies_df.rename(columns={"checkm_completeness":"gtdb_checkm_completeness", "checkm_contamination":"gtdb_checkm_contamination"}, inplace=True)
assemblies_df.drop(columns=["ncbi_genbank_assembly_accession"], inplace=True)
assemblies_df.loc[assemblies_df["assembly_accession"] == "GCA_013112645.1", "ncbi_organism_name"] = "Candidatus Scalindua erythraensis"
column_order = [
        'assembly_accession',
        'ncbi_organism_name',
        'ncbi_strain',
        'ncbi_isolate',
        'ncbi_total_genome_size',
        'ncbi_gc_content',
        'ncbi_assembly_level',
        'ncbi_contig_count',
        'ncbi_taxonomy',
        'gtdb_taxonomy',
        'prokka_annotation',
        'rpoB_identifier',
        'rpoB_hmm_evalue',
        'rpoB_all',
        'rpoB_wo_sporious_taxa',
        'rpoB_extended',
        'rpoB_high_quality',
        'outgroup',
        'gtdb_checkm_completeness',
        'gtdb_checkm_contamination'
]
assemblies_df[column_order].to_csv(
    sys.argv[12], index=False, sep="\t"
)
