
results = "../../processed-data/rpoB-phylogeny/"
data = "../../data/rpoB-phylogeny/"
backup_dir = "/media/argos-emiha442/emiha442/1_projects/1_221123_anammox_pathway_evo/processed-data/rpoB-phylogeny/"
envs = "../envs/"

rule all:
    input:
        #"../processed_data/rpoB-phylogeny/anammox-all-rpoB/anammox-all-rpoB.trimal.initial.treefile",
        expand(results + "anammox-all-rpoB/anammox-pvc-outgroup-rpoB.{cleaning}.trimal.treefile", cleaning=["wo-sporious-taxa", "wo-long-branches"]),
        results + "anammox-extended-rpoB/anammox-extended-outgroup-rpoB.trimal.treefile",
        results + "anammox-hq-rpoB/anammox-hq-outgroup-rpoB.trimal.treefile",
        results + "rpoB-tree-annotation.tsv",

        # Comment out following lines if backup not needed
        #backup_dir + "anammox-all-rpoB/anammox-all-rpoB.trimal.initial.treefile",
        expand(backup_dir + "anammox-all-rpoB/anammox-pvc-outgroup-rpoB.{cleaning}.trimal.treefile", cleaning=["wo-sporious-taxa", "wo-long-branches"]),
        backup_dir + "anammox-extended-rpoB/anammox-extended-outgroup-rpoB.trimal.treefile",
        backup_dir + "anammox-hq-rpoB/anammox-hq-outgroup-rpoB.trimal.treefile",
        backup_dir + "rpoB-tree-annotation.tsv",
        backup_dir + "supplementary-table-1-original.tsv"

#### Create initial rpoB phylogeny ####
rule download_all_anammox_metadata:
    """Download metadata for all assemblies under Candidatus Brocadiia in NCBI taxonomy as of November 2022"""
    output:
        results + "anammox-ncbi-taxonomy.tsv"
    conda:
        envs + "ncbi-datasets.yaml"
    shell:
         """
         datasets summary genome taxon 2517206 --assembly-source GenBank --released-before 11/30/2022 --as-json-lines |\
         dataformat tsv genome --fields accession > {output}
         """

rule download_gtdb_metadata:
    """Download metadata from GTDB"""
    output:
        results + "bac120_metadata_r207.tsv"
    params:
        output_dir=results,
        zip_file=results + "bac120_metadata_r207.tar.gz",
        url="https://data.ace.uq.edu.au/public/gtdb/data/releases/release207/207.0/bac120_metadata_r207.tar.gz"
    shell:
        """
        wget -P {params.output_dir} {params.url} && \
        tar -xvf {params.zip_file} -C {params.output_dir}
        """

rule extract_gtdb_brocadia:
    """Extract all taxa under Brocadiae in the GTDB taxonomy"""
    input:
        results + "bac120_metadata_r207.tsv"
    output:
        results + "anammox-gtdb-taxonomy.tsv"
    shell:
        "python extract-brocadia-from-gtdb.py {input} {output}"

rule merge_ncbi_gtdb_anammox_accessions:
    """Merge the taxa present under Candidatus Brocadiia in NCBI and under Brocadiae in GTDB"""
    input:
        ncbi_taxonomy=results + "anammox-ncbi-taxonomy.tsv",
        gtdb_taxonomy=results + "anammox-gtdb-taxonomy.tsv",
    output:
        results + "anammox-all-accessions.txt",
    shell:
        "python merge-ncbi-gtdb.py {input.ncbi_taxonomy} {input.gtdb_taxonomy} {output}"

rule merge_pvc_anammox_accessions:
    """Merge the anammox assembly accessions with the PVC accessions used as outgroups."""
    input:
        anammox=results + "anammox-all-accessions.txt",
        pvc=data + "pvc-accessions.txt"
    output:
        results + "anammox-pvc-accessions.txt"
    shell:
        "cat {input.anammox} {input.pvc} > {output}"

rule download_assembly_data:
    """ Download genomes and proteomes for all taxa from NCBI GenBank """
    input:
        results + "anammox-pvc-accessions.txt"
    output:
        results + "assembly-data.zip"
    conda:
        envs + "ncbi-datasets.yaml"
    shell:
        "datasets download genome accession --inputfile {input} --include genome,protein --filename {output}"

checkpoint unzip_assembly_data:
    """ Unzip the data """
    input:
        results + "assembly-data.zip"
    output:
        directory(results + "assembly-data")
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
    return expand(os.path.join(results, "prokka", "{assembly_accession}", "{full_accession}.faa"), zip, assembly_accession=assembly_accessions, full_accession=full_accessions)

def ncbi_proteomes(wildcards):
    """ Function to find assemblies which have annotation """
    ck_output = checkpoints.unzip_assembly_data.get(**wildcards).output[0]
    data_folder = os.path.join(ck_output, "ncbi_dataset", "data")
    assembly_accessions = []
    for p in os.listdir(data_folder):
        if os.path.isdir(os.path.join(data_folder, p)):
            if "protein.faa" in os.listdir(os.path.join(data_folder, p)):
                assembly_accessions.append(p)
    return expand(os.path.join(results, "assembly-data/ncbi_dataset/data/", "{assembly_accession}", "protein.faa"), zip, assembly_accession=assembly_accessions)

rule prokka:
    """ Run Prokka on assemblies without annotation """
    input:
        fasta=results + "assembly-data/ncbi_dataset/data/{assembly_accession}/{full_accession}.fna"
    output:
        results + "prokka/{assembly_accession}/{full_accession}.faa"
    params:
        prefix="{full_accession}",
        outdir=results + "prokka/{assembly_accession}"
    threads:
        32
    conda:
        "../envs/prokka.yaml"
    shell:
        """
        prokka --outdir {params.outdir} --prefix {params.prefix} --cpus {threads} {input} --force
        """

rule merge_proteomes:
    """ Merge the proteomes that was downloaded from NCBI with the proteomes from annotated by Prokka """
    input:
        ncbi=ncbi_proteomes,
        prokka=genomes_to_annotate
    output:
        results + "proteomes.faa"
    shell:
        "cat {input.ncbi} {input.prokka} > {output}"

rule download_hmms:
     """ Download the bacteria HMM-file for RpoB """
     output:
         results + "hmm-models/rpoB.hmm.gz"
     shell:
         "wget -O {output} http://eggnogapi5.embl.de/nog_data/file/hmm/COG0085"

rule extract_hmms:
    """ Unzip the HMM """
    input:
        results + "hmm-models/rpoB.hmm.gz"
    output:
        results + "hmm-models/rpoB.hmm"
    shell:
        "gunzip {input}"

rule search_rpoB:
    """ Use hmmsearch to search for RpoB proteins """
    input:
        seq_file=results + "proteomes.faa",
        hmm_profile=results + "hmm-models/rpoB.hmm"
    output:
        table = results + "all-assemblies-hmm-rpoB.tsv",
    conda:
        envs + "hmmer.yaml"
    threads:
        12
    shell:
        "hmmsearch --tblout {output.table} {input.hmm_profile} {input.seq_file}"

rule create_protein_mappings_file:
    input:
        ncbi=ncbi_proteomes,
        prokka=genomes_to_annotate
    output:
        results + "protein-mappings.tsv"
    conda:
        envs + "biopython.yaml"
    shell:
        "python protein-genome-mapping.py {input.ncbi} {input.prokka} {output}"

rule extract_rpoB:
    """ Extract the best hit to RpoB for each taxa, hit must have e-value less than 1e-20. """
    input:
        protein_mappings = results + "protein-mappings.tsv",
        all_proteins = results + "proteomes.faa",
        hmm_search = results + "all-assemblies-hmm-rpoB.tsv",
    output:
        protein = results + "anammox-all-rpoB/all-assemblies-rpoB.faa",
        summary = results + "all-assemblies-rpoB-summary.tsv"
    conda:
        envs + "biopython.yaml"
    shell:
        """
        python extract-rpoB.py \
            {input.protein_mappings} \
            {input.all_proteins} \
            {input.hmm_search} \
            {output.protein} \
            {output.summary}
        """

rule download_ecoli_rpob:
    output:
        results + "outgroup-rpoB/ecoli.faa"
    shell:
        "wget -O {output} https://rest.uniprot.org/uniprotkb/P0A8V2.fasta"

rule download_bsubtilis_rpob:
    output:
        results + "outgroup-rpoB/bsubtilis.faa"
    shell:
        "wget -O {output} https://rest.uniprot.org/uniprotkb/P37870.fasta"

rule merge_pvc_and_outgroup:
    input:
        pvc = results + "anammox-all-rpoB/all-assemblies-rpoB.faa",
        ecoli = results + "outgroup-rpoB/ecoli.faa",
        bsubtilis = results + "outgroup-rpoB/bsubtilis.faa"
    output:
        results + "anammox-all-rpoB/anammox-pvc-outgroup-rpoB.faa"
    shell:
        "cat {input.pvc} {input.ecoli} {input.bsubtilis} > {output}"

rule all_rpoB_alignment:
    input:
        results + "anammox-all-rpoB/anammox-pvc-outgroup-rpoB.faa"
    output:
        results + "anammox-all-rpoB/anammox-pvc-outgroup-rpoB.aln"
    conda:
        envs + "mafft.yaml"
    threads:
        12
    shell:
        "mafft-linsi --reorder --thread {threads} {input} > {output}"

rule all_rpoB_trim_alignment:
    input:
        results + "anammox-all-rpoB/anammox-pvc-outgroup-rpoB.aln"
    output:
        results + "anammox-all-rpoB/anammox-pvc-outgroup-rpoB.trimal.aln"
    conda:
        envs + "trimal.yaml"
    shell:
        "trimal -automated1 -in {input} -out {output}"

rule all_rpoB_phylogeny:
    input:
        results + "anammox-all-rpoB/anammox-pvc-outgroup-rpoB.trimal.aln"
    output:
        results + "anammox-all-rpoB/anammox-pvc-ougroup-rpoB.trimal.initial.treefile"
    params:
        pre=results + "anammox-all-rpoB/anammox-pvc-outgroup-rpoB.trimal.initial"
    conda:
        envs + "iqtree.yaml"
    threads:
        12
    shell:
        "iqtree2 -s {input} -m LG+G4+F -B 1000 -alrt 1000 -pre {params.pre} -nt {threads} -redo"

### Clean up long branches and redo the phylogeny
rule remove_long_branches:
    """ Based on the phylogeny above, taxa with long branches were identified and removed in this step """
    input:
        seq = results + "anammox-all-rpoB/anammox-pvc-outgroup-rpoB.faa",
        seq_to_remove = data + "long-branches.txt"
    output:
        results + "anammox-all-rpoB/anammox-pvc-outgroup-rpoB.wo-long-branches.faa"
    conda:
        envs + "biopython.yaml"
    shell:
        "python ../general/remove-taxa.py {input.seq} {input.seq_to_remove} {output}"

rule clean_all_spurious_taxa:
    input:
        seq = results + "anammox-all-rpoB/anammox-pvc-outgroup-rpoB.faa",
        seq_to_remove = data + "sporious-taxa.txt"
    output:
        results + "anammox-all-rpoB/anammox-pvc-outgroup-rpoB.wo-sporious-taxa.faa"
    conda:
        envs + "biopython.yaml"
    shell:
        "python ../general/remove-taxa.py {input.seq} {input.seq_to_remove} {output}"

rule cleaned_rpoB_alignment:
    """ Align the RpoB sequences long branching taxa removed """
    input:
        results + "anammox-all-rpoB/anammox-pvc-outgroup-rpoB.{cleaning}.faa"
    output:
        results + "anammox-all-rpoB/anammox-pvc-outgroup-rpoB.{cleaning}.aln"
    conda:
        envs + "mafft.yaml"
    threads:
        12
    shell:
        "mafft-linsi --reorder --thread {threads} {input} > {output}"

rule cleaned_rpoB_trim_alignment:
    """ Trim the alignment """
    input:
        results + "anammox-all-rpoB/anammox-pvc-outgroup-rpoB.{cleaning}.aln"
    output:
        results + "anammox-all-rpoB/anammox-pvc-outgroup-rpoB.{cleaning}.trimal.aln"
    conda:
        envs + "trimal.yaml"
    shell:
        "trimal -automated1 -in {input} -out {output}"

rule cleaned_rpoB_phylogeny:
    """ Calculate a new phylogeny based on the cleaned RpoB sequences """
    input:
        results + "anammox-all-rpoB/anammox-pvc-outgroup-rpoB.{cleaning}.trimal.aln"
    output:
        results + "anammox-all-rpoB/anammox-pvc-outgroup-rpoB.{cleaning}.trimal.treefile"
    params:
        pre = results + "anammox-all-rpoB/anammox-pvc-outgroup-rpoB.{cleaning}.trimal"
    conda:
        envs + "iqtree.yaml"
    threads:
        12
    shell:
        "iqtree2 -s {input} -m LG+G4+F -B 1000 -alrt 1000 -pre {params.pre} -nt {threads} -redo"

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
        extended_accessions = data + "anammox-extended-dataset.txt",
        pvc_accessions = data + "pvc-accessions.txt",
        outgroup_accessions = data + "outgroup-accessions.txt",
        proteins = results + "anammox-all-rpoB/anammox-pvc-outgroup-rpoB.faa"
    output:
        results + "anammox-extended-rpoB/anammox-extended-rpoB.faa"
    conda:
        envs + "biopython.yaml"
    shell:
        "python extract-subset.py {input.extended_accessions} {input.pvc_accessions} {input.outgroup_accessions} {input.proteins} {output}"

rule add_ecoli_bsubtilis_to_extended:
    input:
        pvc = results + "anammox-extended-rpoB/anammox-extended-rpoB.faa",
        ecoli = results + "outgroup-rpoB/ecoli.faa",
        bsubtilis = results + "outgroup-rpoB/bsubtilis.faa"
    output:
        results + "anammox-extended-rpoB/anammox-extended-outgroup-rpoB.faa"
    shell:
        "cat {input.pvc} {input.ecoli} {input.bsubtilis} > {output}"

rule extended_dataset_align:
    input:
        results + "anammox-extended-rpoB/anammox-extended-outgroup-rpoB.faa"
    output:
        results + "anammox-extended-rpoB/anammox-extended-outgroup-rpoB.aln"
    conda:
        envs + "mafft.yaml"
    threads:
        12
    shell:
        "mafft-linsi --reorder --thread {threads} {input} > {output}"

rule extended_dataset_trim:
    """ Trim the alignment """
    input:
        results + "anammox-extended-rpoB/anammox-extended-outgroup-rpoB.aln"
    output:
        results + "anammox-extended-rpoB/anammox-extended-outgroup-rpoB.trimal.aln"
    conda:
        envs + "trimal.yaml"
    shell:
        "trimal -automated1 -in {input} -out {output}"

rule extended_dataset_phylogeny:
    """ Calculate a new phylogeny based on the cleaned RpoB sequences """
    input:
        results + "anammox-extended-rpoB/anammox-extended-outgroup-rpoB.trimal.aln"
    output:
        results + "anammox-extended-rpoB/anammox-extended-outgroup-rpoB.trimal.treefile"
    params:
        pre = results + "anammox-extended-rpoB/anammox-extended-outgroup-rpoB.trimal"
    conda:
        envs + "iqtree.yaml"
    threads:
        12
    shell:
        "iqtree2 -s {input} -m LG+G4+F -B 1000 -alrt 1000 -pre {params.pre} -nt {threads} -redo"

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
        extended_accessions = data + "anammox-hq-dataset.txt",
        pvc_accessions = data + "pvc-subset-accessions.txt",
        outgroup_accessions = data + "outgroup-accessions.txt",
        proteins = results + "anammox-all-rpoB/anammox-pvc-outgroup-rpoB.faa"
    output:
        results + "anammox-hq-rpoB/anammox-hq-rpoB.faa"
    conda:
        envs + "biopython.yaml"
    shell:
        "python extract-subset.py {input.extended_accessions} {input.pvc_accessions} {input.outgroup_accessions} {input.proteins} {output}"

rule add_ecoli_bsubtilis_to_hq:
    input:
        pvc = results + "anammox-hq-rpoB/anammox-hq-rpoB.faa",
        ecoli = results + "outgroup-rpoB/ecoli.faa",
        bsubtilis = results + "outgroup-rpoB/bsubtilis.faa"
    output:
        results + "anammox-hq-rpoB/anammox-hq-outgroup-rpoB.faa"
    shell:
        "cat {input.pvc} {input.ecoli} {input.bsubtilis} > {output}"

rule high_quality_dataset_align:
    input:
        results + "anammox-hq-rpoB/anammox-hq-outgroup-rpoB.faa"
    output:
        results + "anammox-hq-rpoB/anammox-hq-outgroup-rpoB.aln"
    conda:
        envs + "mafft.yaml"
    threads:
        12
    shell:
        "mafft-linsi --reorder --thread {threads} {input} > {output}"

rule high_quality_dataset_trim:
    """ Trim the alignment """
    input:
        results + "anammox-hq-rpoB/anammox-hq-outgroup-rpoB.aln"
    output:
        results + "anammox-hq-rpoB/anammox-hq-outgroup-rpoB.trimal.aln"
    conda:
        envs + "trimal.yaml"
    shell:
        "trimal -automated1 -in {input} -out {output}"

rule high_quality_dataset_phylogeny:
    """ Calculate a new phylogeny based on the cleaned RpoB sequences """
    input:
        results + "anammox-hq-rpoB/anammox-hq-outgroup-rpoB.trimal.aln"
    output:
        results + "anammox-hq-rpoB/anammox-hq-outgroup-rpoB.trimal.treefile"
    params:
        pre = results + "anammox-hq-rpoB/anammox-hq-outgroup-rpoB.trimal"
    conda:
        envs + "iqtree.yaml"
    threads:
        12
    shell:
        "iqtree2 -s {input} -m LG+G4+F -B 1000 -alrt 1000 -pre {params.pre} -nt {threads} -redo"

# Create annotation files for the phylogenies that can be opened in fig trees
rule download_tree_annotation_data:
    input:
        anammox = results + "anammox-all-accessions.txt",
        pvc = data + "pvc-accessions.txt"
    output:
        results + "all-assemblies-metadata.tsv"
    conda:
        envs + "ncbi-datasets.yaml"
    shell:
        """
        cat {input.anammox} {input.pvc} | \
        datasets summary genome accession --inputfile - \
            --assembly-source GenBank \
            --released-before 11/30/2022 \
            --as-json-lines | \
        dataformat tsv genome --fields accession,organism-name,organism-infraspecific-strain,organism-infraspecific-isolate,assmstats-total-sequence-len,assmstats-gc-percent,assminfo-level,assmstats-number-of-contigs > {output}
        """

rule fix_tree_annotation_data:
    input:
        results + "all-assemblies-metadata.tsv"
    output:
        results + "rpoB-tree-annotation.tsv"
    shell:
        "python fix-tree-annotation.py {input} {output}"

rule create_supplementary_table_1:
    input:
        metadata = results + "all-assemblies-metadata.tsv",
        ncbi = results + "anammox-ncbi-taxonomy.tsv",
        gtdb = results + "anammox-gtdb-taxonomy.tsv",
        prokka_dir = results + "prokka/",
        rpoB_summary = results + "all-assemblies-rpoB-summary.tsv",
        pvc_accessions = data + "pvc-accessions.txt",
        extended_dataset = data + "anammox-extended-dataset.txt",
        hq_dataset = data + "anammox-hq-dataset.txt",
        pvc_subset = data + "pvc-subset-accessions.txt",
        rpoB_all = results + "anammox-all-rpoB/anammox-pvc-outgroup-rpoB.wo-long-branches.trimal.aln",
        rpoB_cleaned = results + "anammox-all-rpoB/anammox-pvc-outgroup-rpoB.wo-sporious-taxa.trimal.aln"
    output:
        results + "supplementary-table-1-original.tsv"
    conda:
        envs + "biopython.yaml"
    shell:
        """
        python create-tables.py {input.metadata} \
            {input.ncbi} \
            {input.gtdb} \
            {input.prokka_dir} \
            {input.rpoB_summary} \
            {input.pvc_accessions} \
            {input.extended_dataset} \
            {input.hq_dataset} \
            {input.pvc_subset} \
            {input.rpoB_all} \
            {input.rpoB_cleaned} \
            {output}
        """

# Backup results
rule backup_tree_annotation:
    input:
        results + "rpoB-tree-annotation.tsv"
    output:
        backup_dir + "rpoB-tree-annotation.tsv"
    shell:
        "cp {input} {output}"

rule backup_table:
    input:
        results + "supplementary-table-1-original.tsv"
    output:
        backup_dir + "supplementary-table-1-original.tsv"
    shell:
        "cp {input} {output}"

rule backup_all_dataset_initial_phylogeny:
    input:
        results + "anammox-all-rpoB/anammox-pvc-outgroup-rpoB.trimal.initial.treefile"
    output:
        backup_dir + "anammox-all-rpoB/anammox-pvc-outgroup-rpoB.trimal.initial.treefile"
    shell:
        "cp {input} {output}"

rule backup_all_dataset_cleaned_phylogeny:
    input:
        results + "anammox-all-rpoB/anammox-pvc-outgroup-rpoB.{cleaning}.trimal.treefile"
    output:
        backup_dir + "anammox-all-rpoB/anammox-pvc-outgroup-rpoB.{cleaning}.trimal.treefile"
    shell:
        "cp {input} {output}"

rule backup_hq_dataset_phylogeny:
    input:
        results + "anammox-hq-rpoB/anammox-hq-outgroup-rpoB.trimal.treefile"
    output:
        backup_dir + "anammox-hq-rpoB/anammox-hq-outgroup-rpoB.trimal.treefile"
    shell:
        "cp {input} {output}"

rule backup_extended_dataset_phylogeny:
    input:
        results + "anammox-extended-rpoB/anammox-extended-outgroup-rpoB.trimal.treefile"
    output:
        backup_dir + "anammox-extended-rpoB/anammox-extended-outgroup-rpoB.trimal.treefile"
    shell:
        "cp {input} {output}"
