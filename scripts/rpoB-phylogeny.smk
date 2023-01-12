import os

# Backup results
backup_dir = "/media/argos-emiha442/emiha442/1_projects/1_221123_anammox_pathway_evo/processed_data/rpoB-phylogeny/"

rule all:
    input:
        "../processed_data/rpoB-phylogeny/anammox-all-rpoB/anammox-all-rpoB.trimal.initial.treefile",
        "../processed_data/rpoB-phylogeny/anammox-all-rpoB/anammox-all-rpoB.cleaned.trimal.treefile",
        "../processed_data/rpoB-phylogeny/anammox-extended-rpoB/anammox-extended-rpoB.trimal.treefile",
        "../processed_data/rpoB-phylogeny/anammox-hq-rpoB/anammox-hq-rpoB.trimal.treefile",
        "../processed_data/rpoB-phylogeny/rpoB_tree_annnotation.tsv",

        # Comment out following lines if backup not needed
        backup_dir + "anammox-all-rpoB/anammox-all-rpoB.trimal.initial.treefile",
        backup_dir + "anammox-all-rpoB/anammox-all-rpoB.cleaned.trimal.treefile",
        backup_dir + "anammox-extended-rpoB/anammox-extended-rpoB.trimal.treefile",
        backup_dir + "anammox-hq-rpoB/anammox-hq-rpoB.trimal.treefile",
        backup_dir + "rpoB_tree_annnotation.tsv"

rule download_planctomycetes_metadata:
    """Download metadata for assemblies used as outgroup"""
    input:
        "../data/pvc_accessions.txt"
    output:
        "../processed_data/rpoB-phylogeny/planctomycetes/planctomycetes.metadata.tsv"
    conda:
        "envs/ncbi-datasets.yaml"
    shell:
        """
        datasets summary genome accession --inputfile {input} --as-json-lines |\
        dataformat tsv genome --fields accession,organism-name,assmstats-total-sequence-len,assmstats-number-of-contigs,assminfo-biosample-attribute-name,assminfo-biosample-attribute-value > {output}
        """

rule download_anammox_metadata:
    input:
        "../data/anammox_accessions.txt"
    output:
        "../processed_data/rpoB-phylogeny/anammox.metadata.tsv"
    conda:
        "envs/ncbi-datasets.yaml"
    shell:
        """
        datasets summary genome accession --inputfile {input} --as-json-lines |\
        dataformat tsv genome --fields accession,organism-name,assmstats-total-sequence-len,assmstats-number-of-contigs,assminfo-biosample-attribute-name,assminfo-biosample-attribute-value > {output}
        """
#### Create initial rpoB phylogeny ####
rule download_all_anammox_metadata:
    """Download metadata for all assemblies under Candidatus Brocadiia in NCBI taxonomy as of November 2022"""
    output:
        "../processed_data/rpoB-phylogeny/anammox.ncbi_taxonomy.tsv"
    conda:
        "envs/ncbi-datasets.yaml"
    shell:
         """
         datasets summary genome taxon 2517206 --assembly-source GenBank --released-before 30/11/2022 --as-json-lines |\
         dataformat tsv genome --fields accession,organism-name,assmstats-total-sequence-len,assmstats-number-of-contigs,assminfo-biosample-attribute-name,assminfo-biosample-attribute-value > {output}
         """

rule download_gtdb_metadata:
    """Download metadata from GTDB"""
    output:
        "../processed_data/rpoB-phylogeny/bac120_metadata_r207.tsv"
    params:
        output_dir="../processed_data/rpoB-phylogeny/",
        zip_file="../processed_data/rpoB-phylogeny/bac120_metadata_r207.tar.gz",
        url="https://data.ace.uq.edu.au/public/gtdb/data/releases/release207/207.0/bac120_metadata_r207.tar.gz"
    shell:
        """
        wget -P {params.output_dir} {params.url} && \
        tar -xvf {params.zip_file} -C {params.output_dir}
        """

rule extract_gtdb_brocadia:
    """Extract all taxa under Brocadiae in the GTDB taxonomy"""
    input:
        "../processed_data/rpoB-phylogeny/bac120_metadata_r207.tsv"
    output:
        "../processed_data/rpoB-phylogeny/anammox.gtdb_taxonomy.tsv"
    shell:
        "python rpoB-phylogeny-gtdb-brocadia.py {input} {output}"

rule merge_ncbi_gtdb_anammox_accessions:
    """Merge the taxa present under Candidatus Brocadiia in NCBI and under Brocadiae in GTDB"""
    input:
        ncbi_taxonomy="../processed_data/rpoB-phylogeny/anammox.ncbi_taxonomy.tsv",
        gtdb_taxonomy="../processed_data/rpoB-phylogeny/anammox.gtdb_taxonomy.tsv",
    output:
        "../processed_data/rpoB-phylogeny/anammox.all_accessions.txt",
    shell:
        "python rpoB-phylogeny-merge-ncbi-gtdb.py {input.ncbi_taxonomy} {input.gtdb_taxonomy} {output}"

rule merge_pvc_anammox_accessions:
    """Merge the anammox assembly accessions with the PVC accessions used as outgroups."""
    input:
        anammox="../processed_data/rpoB-phylogeny/anammox.all_accessions.txt",
        pvc="../data/pvc_accessions.txt"
    output:
        "../processed_data/rpoB-phylogeny/anammox_pvc_accessions.txt"
    shell:
        "cat {input.anammox} {input.pvc} > {output}"

rule download_assembly_data:
    """ Download genomes and proteomes for all taxa from NCBI GenBank """
    input:
        "../processed_data/rpoB-phylogeny/anammox_pvc_accessions.txt"
    output:
        "../processed_data/rpoB-phylogeny/assembly_data.zip"
    conda:
        "envs/ncbi-datasets.yaml"
    shell:
        "datasets download genome accession --inputfile {input} --include genome,protein --filename {output}"

checkpoint unzip_assembly_data:
    """ Unzip the data """
    input:
        "../processed_data/rpoB-phylogeny/assembly_data.zip"
    output:
        directory("../processed_data/rpoB-phylogeny/assembly_data")
    shell:
        "unzip {input} -d {output}"

def genomes_to_annotate(wildcards):
    """ Function to find assemblies without annotation in the downloaded dataset """
    ck_output = checkpoints.unzip_assembly_data.get(**wildcards).output[0]
    data_folder = os.path.join(ck_output, "ncbi_dataset", "data")
    full_accessions = []
    assembly_accessions = []
    for p in os.listdir(data_folder):
        if os.path.isdir(os.path.join(data_folder, p)):
            if "protein.faa" not in os.listdir(os.path.join(data_folder, p)):
                assembly_accessions.append(p)
                full_accession = os.listdir(os.path.join(data_folder, p))[0]
                full_accession = os.path.splitext(full_accession)[0]
                full_accessions.append(full_accession)
    return expand(os.path.join("../processed_data/rpoB-phylogeny/prokka/", "{assembly_accession}", "{full_accession}.faa"), zip, assembly_accession=assembly_accessions, full_accession=full_accessions)

def ncbi_proteomes(wildcards):
    """ Function to find assemblies which have annotation """
    ck_output = checkpoints.unzip_assembly_data.get(**wildcards).output[0]
    data_folder = os.path.join(ck_output, "ncbi_dataset", "data")
    assembly_accessions = []
    for p in os.listdir(data_folder):
        if os.path.isdir(os.path.join(data_folder, p)):
            if "protein.faa" in os.listdir(os.path.join(data_folder, p)):
                assembly_accessions.append(p)
    return expand(os.path.join("../processed_data/rpoB-phylogeny/assembly_data/ncbi_dataset/data/", "{assembly_accession}", "protein.faa"), zip, assembly_accession=assembly_accessions)

rule prokka:
    """ Run Prokka on assemblies without annotation """
    input:
        fasta="../processed_data/rpoB-phylogeny/assembly_data/ncbi_dataset/data/{assembly_accession}/{full_accession}.fna"
    output:
        "../processed_data/rpoB-phylogeny/prokka/{assembly_accession}/{full_accession}.faa"
    params:
        prefix="{full_accession}",
        outdir="../processed_data/rpoB-phylogeny/prokka/{assembly_accession}"
    threads:
        24
    conda:
        "envs/prokka.yaml"
    shell:
        """
        prokka --outdir {params.outdir} --prefix {params.prefix} --cpus {threads} {input} --kingdom bacteria --force
        """

rule merge_proteomes:
    """ Merge the proteomes that was downloaded from NCBI with the proteomes from annotated by Prokka """
    input:
        ncbi=ncbi_proteomes,
        prokka=genomes_to_annotate
    output:
        temp("../processed_data/rpoB-phylogeny/proteomes.faa")
    shell:
        "cat {input.ncbi} {input.prokka} > {output}"

rule download_hmms:
     """ Download the bacteria HMM-file for RpoB """
     output:
         "../processed_data/rpoB-phylogeny/hmm-models/rpoB.hmm.gz"
     shell:
         "wget -O {output} http://eggnogapi5.embl.de/nog_data/file/hmm/COG0085"

rule extract_hmms:
    """ Unzip the HMM """
    input:
        "../processed_data/rpoB-phylogeny/hmm-models/rpoB.hmm.gz"
    output:
        "../processed_data/rpoB-phylogeny/hmm-models/rpoB.hmm"
    shell:
        "gunzip {input}"

rule search_rpoB:
    """ Use hmmsearch to search for RpoB proteins """
    input:
        seq_file="../processed_data/rpoB-phylogeny/proteomes.faa",
        hmm_profile="../processed_data/rpoB-phylogeny/hmm-models/rpoB.hmm"
    output:
        table = "../processed_data/rpoB-phylogeny/all_assemblies.rpoB.tsv",
    conda:
        "envs/hmmer.yaml"
    threads:
        12
    shell:
        "hmmsearch --tblout {output.table} {input.hmm_profile} {input.seq_file}"

rule create_protein_mappings_file:
    input:
        ncbi=ncbi_proteomes,
        prokka=genomes_to_annotate
    output:
        temp("../processed_data/rpoB-phylogeny/protein_mappings.tsv")
    conda:
        "envs/biopython.yaml"
    shell:
        "python rpoB-phylogeny-protein-genome-mapping.py {input.ncbi} {input.prokka} {output}"

rule extract_rpoB:
    """ Extract the best hit to RpoB for each taxa, hit must have e-value less than 1e-20. """
    input:
        protein_mappings="../processed_data/rpoB-phylogeny/protein_mappings.tsv",
        all_proteins="../processed_data/rpoB-phylogeny/proteomes.faa",
        hmm_search="../processed_data/rpoB-phylogeny/all_assemblies.rpoB.tsv",
    output:
        "../processed_data/rpoB-phylogeny/anammox-all-rpoB/anammox-all-rpoB.faa"
    conda:
        "envs/biopython.yaml"
    shell:
        """
        python rpoB-phylogeny-extract-rpoB.py \
            {input.protein_mappings} \
            {input.all_proteins} \
            {input.hmm_search} \
            {output}
        """
rule all_rpoB_alignment:
    input:
        "../processed_data/rpoB-phylogeny/anammox-all-rpoB/anammox-all-rpoB.faa"
    output:
        "../processed_data/rpoB-phylogeny/anammox-all-rpoB/anammox-all-rpoB.aln"
    conda:
        "envs/mafft.yaml"
    threads:
        12
    shell:
        "mafft-linsi --reorder --thread {threads} {input} > {output}"

rule all_rpoB_trim_alignment:
    input:
        "../processed_data/rpoB-phylogeny/anammox-all-rpoB/anammox-all-rpoB.aln"
    output:
        "../processed_data/rpoB-phylogeny/anammox-all-rpoB/anammox-all-rpoB.trimal.aln"
    conda:
        "envs/trimal.yaml"
    shell:
        "trimal -automated1 -in {input} -out {output}"


rule all_rpoB_phylogeny:
    input:
        "../processed_data/rpoB-phylogeny/anammox-all-rpoB/anammox-all-rpoB.trimal.aln"
    output:
        "../processed_data/rpoB-phylogeny/anammox-all-rpoB/anammox-all-rpoB.trimal.initial.treefile"
    params:
        pre="../processed_data/rpoB-phylogeny/anammox-all-rpoB/anammox-all-rpoB.trimal.initial"
    conda:
        "envs/iqtree.yaml"
    threads:
        12
    shell:
        "iqtree2 -s {input} -m LG+G4+F -B 1000 -pre {params.pre} -nt {threads} -redo"

### Clean up long branches and redo the phylogeny
rule remove_long_branches:
    """ Based on the phylogeny above, taxa with long branches were identified and removed in this step """
    input:
        "../processed_data/rpoB-phylogeny/anammox-all-rpoB/anammox-all-rpoB.faa"
    output:
        "../processed_data/rpoB-phylogeny/anammox-all-rpoB/anammox-all-rpoB.cleaned.faa"
    conda:
        "envs/biopython.yaml"
    shell:
        "python rpoB-phylogeny-remove-long-branches.py {input} {output}"

rule cleaned_rpoB_alignment:
    """ Align the RpoB sequences long branching taxa removed """
    input:
        "../processed_data/rpoB-phylogeny/anammox-all-rpoB/anammox-all-rpoB.cleaned.faa"
    output:
        "../processed_data/rpoB-phylogeny/anammox-all-rpoB/anammox-all-rpoB.cleaned.aln"
    conda:
        "envs/mafft.yaml"
    threads:
        12
    shell:
        "mafft-linsi --reorder --thread {threads} {input} > {output}"

rule cleaned_rpoB_trim_alignment:
    """ Trim the alignment """
    input:
        "../processed_data/rpoB-phylogeny/anammox-all-rpoB/anammox-all-rpoB.cleaned.aln"
    output:
        "../processed_data/rpoB-phylogeny/anammox-all-rpoB/anammox-all-rpoB.cleaned.trimal.aln"
    conda:
        "envs/trimal.yaml"
    shell:
        "trimal -automated1 -in {input} -out {output}"

rule cleaned_rpoB_phylogeny:
    """ Calculate a new phylogeny based on the cleaned RpoB sequences """
    input:
        "../processed_data/rpoB-phylogeny/anammox-all-rpoB/anammox-all-rpoB.cleaned.trimal.aln"
    output:
        "../processed_data/rpoB-phylogeny/anammox-all-rpoB/anammox-all-rpoB.cleaned.trimal.treefile"
    params:
        pre="../processed_data/rpoB-phylogeny/anammox-all-rpoB/anammox-all-rpoB.cleaned.trimal"
    conda:
        "envs/iqtree.yaml"
    threads:
        12
    shell:
        "iqtree2 -s {input} -m LG+G4+F -B 1000 -pre {params.pre} -nt {threads} -redo"

# Based on the RpoB phylogeny produced above and on genome assembly information
# in supplementary table 1 a set of high quality genomes were selected for
# in-depth analysis. Also an extended set was produced which includes additional taxa
# including recently described early-diverging anammox clade. Below is the workflow
# to produce RpoB-phylogenies for the two different datasets.

#### RpoB phylogeny of the extended anammox dataset ####
rule extended_dataset_extract_rpoB:
    """
    Accessions for the anammox genomes included in the extended dataset is
    stored in text file in the datafolder. Additional information about this
    taxa can be found in supplementary table 3.
    This step extract the RpoB sequences for the anammox accessions in the extended
    dataset and also the accessions in the PVC dataset.
    """
    input:
        extended_accessions="../data/anammox_extended_dataset.txt",
        pvc_accessions="../data/pvc_accessions.txt",
        proteins="../processed_data/rpoB-phylogeny/anammox-all-rpoB/anammox-all-rpoB.faa"
    output:
        "../processed_data/rpoB-phylogeny/anammox-extended-rpoB/anammox-extended-rpoB.faa"
    conda:
        "envs/biopython.yaml"
    shell:
        "python rpoB-phylogeny-extract-subset.py {input.extended_accessions} {input.pvc_accessions} {input.proteins} {output}"

rule extended_dataset_align:
    input:
        "../processed_data/rpoB-phylogeny/anammox-extended-rpoB/anammox-extended-rpoB.faa"
    output:
        "../processed_data/rpoB-phylogeny/anammox-extended-rpoB/anammox-extended-rpoB.aln"
    conda:
        "envs/mafft.yaml"
    threads:
        12
    shell:
        "mafft-linsi --reorder --thread {threads} {input} > {output}"

rule extended_dataset_trim:
    """ Trim the alignment """
    input:
        "../processed_data/rpoB-phylogeny/anammox-extended-rpoB/anammox-extended-rpoB.aln"
    output:
        "../processed_data/rpoB-phylogeny/anammox-extended-rpoB/anammox-extended-rpoB.trimal.aln"
    conda:
        "envs/trimal.yaml"
    shell:
        "trimal -automated1 -in {input} -out {output}"

rule extended_dataset_phylogeny:
    """ Calculate a new phylogeny based on the cleaned RpoB sequences """
    input:
        "../processed_data/rpoB-phylogeny/anammox-extended-rpoB/anammox-extended-rpoB.trimal.aln"
    output:
        "../processed_data/rpoB-phylogeny/anammox-extended-rpoB/anammox-extended-rpoB.trimal.treefile"
    params:
        pre="../processed_data/rpoB-phylogeny/anammox-extended-rpoB/anammox-extended-rpoB.trimal"
    conda:
        "envs/iqtree.yaml"
    threads:
        12
    shell:
        "iqtree2 -s {input} -m LG+G4+F -B 1000 -pre {params.pre} -nt {threads} -redo"

#### RpoB phylogeny of the high quality anammox dataset ####
rule high_quality_dataset_extract_rpoB:
    """
    Accessions for the anammox genomes included in the high-quality dataset are
    stored in text file in the datafolder. Additional information about these
    taxa can be found in supplementary table 2.
    This step extract the RpoB sequences for the accessions in the high quality dataset
    and accesstions for the PVC dataset.
    """
    input:
        extended_accessions="../data/anammox_hq_dataset.txt",
        pvc_accessions="../data/pvc_accessions.txt",
        proteins="../processed_data/rpoB-phylogeny/anammox-all-rpoB/anammox-all-rpoB.faa"
    output:
        "../processed_data/rpoB-phylogeny/anammox-hq-rpoB/anammox-hq-rpoB.faa"
    conda:
        "envs/biopython.yaml"
    shell:
        "python rpoB-phylogeny-extract-subset.py {input.extended_accessions} {input.pvc_accessions} {input.proteins} {output}"

rule high_quality_dataset_align:
    input:
        "../processed_data/rpoB-phylogeny/anammox-hq-rpoB/anammox-hq-rpoB.faa"
    output:
        "../processed_data/rpoB-phylogeny/anammox-hq-rpoB/anammox-hq-rpoB.aln"
    conda:
        "envs/mafft.yaml"
    threads:
        12
    shell:
        "mafft-linsi --reorder --thread {threads} {input} > {output}"

rule high_quality_dataset_trim:
    """ Trim the alignment """
    input:
        "../processed_data/rpoB-phylogeny/anammox-hq-rpoB/anammox-hq-rpoB.aln"
    output:
        "../processed_data/rpoB-phylogeny/anammox-hq-rpoB/anammox-hq-rpoB.trimal.aln"
    conda:
        "envs/trimal.yaml"
    shell:
        "trimal -automated1 -in {input} -out {output}"

rule high_quality_dataset_phylogeny:
    """ Calculate a new phylogeny based on the cleaned RpoB sequences """
    input:
        "../processed_data/rpoB-phylogeny/anammox-hq-rpoB/anammox-hq-rpoB.trimal.aln"
    output:
        "../processed_data/rpoB-phylogeny/anammox-hq-rpoB/anammox-hq-rpoB.trimal.treefile"
    params:
        pre="../processed_data/rpoB-phylogeny/anammox-hq-rpoB/anammox-hq-rpoB.trimal"
    conda:
        "envs/iqtree.yaml"
    threads:
        12
    shell:
        "iqtree2 -s {input} -m LG+G4+F -B 1000 -pre {params.pre} -nt {threads} -redo"

# Create annotation files for the phylogenies that can be opened in fig trees
rule download_tree_annotation_data:
    input:
        anammox="../processed_data/rpoB-phylogeny/anammox.all_accessions.txt",
        pvc="../data/pvc_accessions.txt"
    output:
        "../processed_data/rpoB-phylogeny/tree_annnotation.raw.tsv"
    conda:
        "envs/ncbi-datasets.yaml"
    shell:
        """
        cat {input.anammox} {input.pvc} | \
        datasets summary genome accession --inputfile - \
            --assembly-source GenBank \
            --released-before 11/30/2022 \
            --as-json-lines | \
        dataformat tsv genome --fields accession,organism-name,organism-infraspecific-strain,organism-infraspecific-isolate > {output}
        """
rule fix_tree_annotation_data:
    input:
        "../processed_data/rpoB-phylogeny/tree_annnotation.raw.tsv"
    output:
        "../processed_data/rpoB-phylogeny/rpoB_tree_annnotation.tsv"
    shell:
        "python rpoB-phylogeny-fix-tree-annotation.py {input} {output}"

# Backup results
rule backup_tree_annotation:
    input:
        "../processed_data/rpoB-phylogeny/rpoB_tree_annnotation.tsv"
    output:
        backup_dir + "rpoB_tree_annnotation.tsv"

rule backup_all_dataset_initial_phylogeny:
    input:
        "../processed_data/rpoB-phylogeny/anammox-all-rpoB/anammox-all-rpoB.trimal.initial.treefile"
    output:
        backup_dir + "anammox-all-rpoB/anammox-all-rpoB.trimal.initial.treefile"
    shell:
        "cp {input} {output}"

rule backup_all_dataset_cleaned_phylogeny:
    input:
        "../processed_data/rpoB-phylogeny/anammox-all-rpoB/anammox-all-rpoB.cleaned.trimal.treefile"
    output:
        backup_dir + "anammox-all-rpoB/anammox-all-rpoB.cleaned.trimal.treefile"
    shell:
        "cp {input} {output}"

rule backup_hq_dataset_phylogeny:
    input:
        "../processed_data/rpoB-phylogeny/anammox-hq-rpoB/anammox-hq-rpoB.trimal.treefile"
    output:
        backup_dir + "anammox-hq-rpoB/anammox-hq-rpoB.trimal.treefile"
    shell:
        "cp {input} {output}"

rule backup_extended_dataset_phylogeny:
    input:
        "../processed_data/rpoB-phylogeny/anammox-extended-rpoB/anammox-extended-rpoB.trimal.treefile"
    output:
        backup_dir + "anammox-extended-rpoB/anammox-extended-rpoB.trimal.treefile"
    shell:
        "cp {input} {output}"
