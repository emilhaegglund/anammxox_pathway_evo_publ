backup_dir = "/media/argos-emiha442/emiha442/1_projects/1_221123_anammox_pathway_evo/processed_data/hzs-outside-anammox/"
results = "../processed_data/hzs-outside-anammox/"
SUBUNITS = ["hzs_a", "hzs_bc"]
rule all:
    input:
        results + "hzs.gtdb.w_taxa.tsv.gz",
        expand(results + "{subunit}.gtdb.refseq.hits.trimal.fasttree", subunit=SUBUNITS),
        expand(results + "{subunit}.gtdb.refseq.hits.gap_remove_1.trimal.fasttree", subunit=SUBUNITS),
        expand(results + "{subunit}.gtdb.refseq.hits.remove_long_branches.trimal.fasttree", subunit=SUBUNITS),
        expand(results + "{subunit}.gtdb.refseq.hits.remove_long_branches.trimal.mad_root.treefile", subunit=SUBUNITS),
        expand(results + "{subunit}.blast.structure.genome.heme.annotation.tsv", subunit=SUBUNITS),
        expand(results + "{subunit}.gtdb.hits.faa", subunit=SUBUNITS),
        results + "hzs_bc.sister_group.anammox.aln",
        results + "hzs.blast.summary.svg",
        results + "hzs_bc_refseqs_proteomes.faa",
        results + "hzs_bc_refseqs_proteomes.hzs_a_domanin.filtered.tsv",

rule download_queries:
    "Download HZS-B, HZS-C, and HZS-A from Kuenenia stuttgartiensis and HZS-BC and HZS-A from Scalinduae japonica"
    output:
        results + "query.faa"
    conda:
        "envs/entrez.yaml"
    shell:
        "efetch -id SOH05198.1,SOH05199.1,SOH05200.1,GAX62881.1,GAX62882.1 -db protein -format fasta > {output}"

rule blast_gtdb:
    "Search for homologous proteins among the representative species in GTDB v207"
    input:
        results + "query.faa"
    output:
        results + "hzs.gtdb.tsv.gz"
    params:
        db="../data/gtdb_representatives_database/gtdb_representatives.dmnd"
    conda:
        "envs/diamond.yaml"
    threads:
        32
    shell:
        """
        diamond blastp --db {params.db} \
            --query {input} \
            --outfmt 6 qseqid sseqid pident positive mismatch gapopen qlen slen length qstart qend sstart send evalue bitscore stitle full_sseq \
            --out {output} \
            --very-sensitive \
            --max-target-seqs 5000 \
            --min-score 80 \
            --compress 1 \
            --threads {threads} \
        """

rule add_accessions:
    "Add the accession number of the assembly for each of the identified proteins"
    input:
        blast=results + "hzs.gtdb.tsv.gz",
        prot="../data/gtdb_representatives_database/prot2accession.tsv.gz"
    output:
        results + "hzs.gtdb.w_accessions.tsv.gz"
    shell:
        "python hzs-outside-anammox-merge-blast-accessions.py {input.blast} {input.prot} {output}"

rule add_taxa:
    "Add taxonomic information for each of the identified proteins"
    input:
        blast=results + "hzs.gtdb.w_accessions.tsv.gz",
        gtdb="../data/gtdb_representatives_database/gtdb_representatives.sample_gtdb.metadata.tsv",
    output:
        results + "hzs.gtdb.w_taxa.tsv.gz"
    shell:
        "python hzs-outside-anammox-merge-blast-taxa.py {input.blast} {input.gtdb} {output}"

rule plot_summary_blast_results:
    "Plot a summary of the blast results. Presented in supplementary figure X"
    input:
        results + "hzs.gtdb.w_taxa.tsv.gz"
    output:
        results + "hzs.blast.summary.svg"
    shell:
        "python plot_blast_results_distribution.py {input} {output}"

# Prepare phylogeny from RefSeq homologs
rule extract_hzs_a_refseq_sequences:
    """Extract protein sequence for all hits to HZS-A for entries derived from
    RefSeq"""
    input:
        blast=results + "hzs.gtdb.w_taxa.tsv.gz"
    output:
        results + "hzs_a.gtdb.refseq.no_anammox.hits.faa"
    shell:
        """python extract-blast-hits.py --blast {input.blast} \
            --fasta {output} \
            --queries SOH05200.1,GAX62882.1 \
            --scover 0.5 \
            --qcover 0.5 \
            --source refseq
        """

rule extract_hzs_bc_refseq_sequences:
    """Extract protein sequence for all hits to HZS-BC, HZS-B, and HZS-C for
    entries derived from RefSeq."""
    input:
        blast=results + "hzs.gtdb.w_taxa.tsv.gz"
    output:
        results + "hzs_bc.gtdb.refseq.no_anammox.hits.faa"
    shell:
        """python extract-blast-hits.py --blast {input.blast} \
            --fasta {output} \
            --queries SOH05198.1,SOH05199.1,GAX62881.1 \
            --scover 0.5 \
            --qcover 0.5 \
            --source refseq
        """

rule merge_anammox_and_no_anammox_refseqs:
    """Protein sequence for HZS-A, HZS-B, and HZS-C for anammox sequences from
    the HQ-dataset is located in the data folder. For B and C the sequences have
    been manually fused. Only one copy for each of the species have been
    included. Merge the anammox protein sequences with the protein sequences
    for the RefSeq hits."""
    input:
        no_anammox=results + "{subunit}.gtdb.refseq.no_anammox.hits.faa",
        anammox="../data/hzs-outside-anammox/{subunit}.faa"
    output:
        results + "{subunit}.gtdb.refseq.hits.faa"
    shell:
        "cat {input.anammox} {input.no_anammox} > {output}"

rule align_refseqs:
    "Align the sequences using mafft."
    input:
        results + "{subunit}.gtdb.refseq.hits.faa"
    output:
        results + "{subunit}.gtdb.refseq.hits.aln"
    conda:
        "./envs/mafft.yaml"
    threads:
        12
    shell:
        "mafft-linsi --reorder --thread {threads} {input} > {output}"

rule trim_refseqs:
    "Trim the aligment using TrimAl."
    input:
        results + "{subunit}.gtdb.refseq.hits.aln"
    output:
        results + "{subunit}.gtdb.refseq.hits.trimal.aln"
    conda:
        "./envs/trimal.yaml"
    threads:
        12
    shell:
        "trimal -automated1 -keepheader -in {input} -out {output}"

rule fasttree_refseqs:
    "Run FastTree to get first overview of data."
    input:
        results + "{subunit}.gtdb.refseq.hits.trimal.aln"
    output:
        results + "{subunit}.gtdb.refseq.hits.trimal.fasttree"
    conda:
        "./envs/fasttree.yaml"
    shell:
        "fasttree -lg {input} > {output}"

rule remove_gappy_sequences_refseqs_1:
    input:
        aln=results + "{subunit}.gtdb.refseq.hits.trimal.aln",
        seq=results + "{subunit}.gtdb.refseq.hits.faa"
    output:
        results + "{subunit}.gtdb.refseq.hits.gap_remove_1.faa"
    conda:
        "./envs/biopython.yaml"
    shell:
        "python remove_more_50_gaps.py {input.aln} {input.seq} {output}"

rule align_refseqs_after_gap_remove_1:
    input:
        results + "{subunit}.gtdb.refseq.hits.gap_remove_1.faa"
    output:
        results + "{subunit}.gtdb.refseq.hits.gap_remove_1.aln"
    conda:
        "./envs/mafft.yaml"
    threads:
        12
    shell:
        "mafft-linsi --reorder --thread {threads} {input} > {output}"

rule trim_refseqs_after_gap_remove_1:
    input:
        results + "{subunit}.gtdb.refseq.hits.gap_remove_1.aln"
    output:
        results + "{subunit}.gtdb.refseq.hits.gap_remove_1.trimal.aln"
    conda:
        "./envs/trimal.yaml"
    threads:
        12
    shell:
        "trimal -automated1 -keepheader -in {input} -out {output}"

rule fasttree_refseqs_after_gap_remove_1:
    input:
        results + "{subunit}.gtdb.refseq.hits.gap_remove_1.trimal.aln"
    output:
        results + "{subunit}.gtdb.refseq.hits.gap_remove_1.trimal.fasttree"
    conda:
        "./envs/fasttree.yaml"
    shell:
        "fasttree -lg {input} > {output}"

rule remove_long_branches:
    input:
        sequences=results + "{subunit}.gtdb.refseq.hits.gap_remove_1.faa",
        accessions="../data/hzs-outside-anammox/{subunit}.long_branches.accessions.txt"
    output:
        results + "{subunit}.gtdb.refseq.hits.remove_long_branches.faa"
    conda:
        "./envs/biopython.yaml"
    shell:
        "python rpoB-phylogeny-remove-taxa.py {input.sequences} {input.accessions} {output}"

rule align_refseqs_after_remove_long_branches:
    input:
        results + "{subunit}.gtdb.refseq.hits.remove_long_branches.faa"
    output:
        results + "{subunit}.gtdb.refseq.hits.remove_long_branches.aln"
    conda:
        "./envs/mafft.yaml"
    threads:
        12
    shell:
        "mafft-linsi --reorder --thread {threads} {input} > {output}"

rule trim_refseqs_after_remove_long_branches:
    input:
        results + "{subunit}.gtdb.refseq.hits.remove_long_branches.aln"
    output:
        results + "{subunit}.gtdb.refseq.hits.remove_long_branches.trimal.aln"
    conda:
        "./envs/trimal.yaml"
    threads:
        12
    shell:
        "trimal -automated1 -keepheader -in {input} -out {output}"

rule fasttree_refseqs_after_remove_long_branches:
    input:
        results + "{subunit}.gtdb.refseq.hits.remove_long_branches.trimal.aln"
    output:
        results + "{subunit}.gtdb.refseq.hits.remove_long_branches.trimal.fasttree"
    conda:
        "./envs/fasttree.yaml"
    shell:
        "fasttree -lg {input} > {output}"

rule iqtree_refseqs:
    input:
        results + "{subunit}.gtdb.refseq.hits.remove_long_branches.trimal.aln"
    output:
        results + "{subunit}.gtdb.refseq.hits.remove_long_branches.trimal.treefile"
    conda:
        "envs/iqtree.yaml"
    params:
        pre=results + "{subunit}.gtdb.refseq.hits.remove_long_branches.trimal"
    threads:
        12
    shell:
        "iqtree2 -s {input} -mset LG,WAG -alrt 1000 -B 1000 -nt {threads} -pre {params.pre} -redo"

rule mad_root:
    input:
        results + "{subunit}.gtdb.refseq.hits.remove_long_branches.trimal.treefile"
    output:
        results + "{subunit}.gtdb.refseq.hits.remove_long_branches.trimal.mad_root.treefile"
    shell:
        "../bin/MADroot/madRoot {input} > {output}"

# Align only Brocadiae and sister group
rule extract_sister_group:
    input:
        accessions="../data/hzs-outside-anammox/hzs_bc.sister_group.accessions.txt",
        sequences=results + "hzs_bc.gtdb.refseq.hits.faa"
    output:
        results + "hzs_bc.sister_group.faa"
    conda:
        "./envs/seqtk.yaml"
    shell:
        "seqtk subseq {input.sequences} {input.accessions} > {output}"

rule merge_sister_group_anammox:
    input:
        sister_group=results + "hzs_bc.sister_group.faa",
        anammox="../data/hzs-outside-anammox/hzs_bc.faa"
    output:
        results + "hzs_bc.sister_group.anammox.faa"
    shell:
        "cat {input.sister_group} {input.anammox} > {output}"

rule align_anammox_sister_group:
    input:
        results + "hzs_bc.sister_group.anammox.faa"
    output:
        results + "hzs_bc.sister_group.anammox.aln"
    conda:
        "./envs/mafft.yaml"
    shell:
        "mafft-linsi --reorder {input} > {output}"

# Phylogeny based on clustering of refseqs
rule cluster_refseqs:
    input:
        results + "/{subunit}.gtdb.refseq.hits.gap_remove_1.faa"
    output:
        results + "/{subunit}.gtdb.refseq.hits.gap_remove_1.clustered.faa"
    conda:
        "./envs/cd-hit.yaml"
    shell:
        "cd-hit -i {input} -o {output} -c 0.9"

rule align_cluster_refseqs:
    input:
        results + "/{subunit}.gtdb.refseq.hits.gap_remove_1.clustered.faa"
    output:
        results + "/{subunit}.gtdb.refseq.hits.gap_remove_1.clustered.aln"
    conda:
        "./envs/mafft.yaml"
    threads:
        12
    shell:
        "mafft-linsi --reorder --thread {threads} {input} > {output}"

rule trim_clustered_refseqs:
    input:
        results + "/{subunit}.gtdb.refseq.hits.gap_remove_1.clustered.aln"
    output:
        results + "/{subunit}.gtdb.refseq.hits.gap_remove_1.clustered.trimal.aln"
    conda:
        "./envs/trimal.yaml"
    shell:
        "trimal -automated1 -in {input} -out {output}"

rule iqtree_clustered_refseqs:
    input:
        results + "/{subunit}.gtdb.refseq.hits.gap_remove_1.clustered.trimal.aln"
    output:
        results + "/{subunit}.gtdb.refseq.hits.gap_remove_1.clustered.trimal.treefile"
    conda:
        "envs/iqtree.yaml"
    params:
        pre=results + "/{subunit}.gtdb.refseq.hits.gap_remove_1.clustered.trimal"
    threads:
        12
    shell:
        "iqtree2 -s {input} -mset LG -bb 1000 -nt {threads} -pre {params.pre} -redo"

# Download GBFF for comparative analysis of the Refseqs. The following steps
# download proteomes and GBK files for all proteomes in RefSeq with a hit to the
# HZS-BC subunit. Since HZS-A is less conserved compared to HZS-BC we search the
# proteomes for the HZS-A middle domain using hmmsearch. The GBKs is also used
# to plot the genome structure around the HZS-BC homologs.
rule extract_refseq_accessions:
    input:
        results + "{subunit}.blast.structure.genome.tsv"
    output:
        results + "{subunit}.gtdb.refseq.hits.accessions.txt"
    params:
        subunit="{subunit}"
    shell:
        "python extract_refseq_accession.py {input} {params.subunit} > {output}"

# The following was run outside the Snakemake workflow.
#rule download_refseq_gbff:
#    input:
#        results + "{subunit}.gtdb.refseq.hits.accessions.txt"
#    output:
#        results + "{subunit}.refseq.gbff.zip"
#    conda:
#        "./envs/ncbi-datasets.yaml"
#    shell:
#        "datasets download genome accession --inputfile {input} --filename {output} --include gbff,protein"

rule unzip_refseq_gbff:
    input:
        results + "{subunit}.refseq.gbff.zip"
    output:
        directory(results + "{subunit}_refseqs_gbff")
    shell:
        "unzip {input} -d {output}"

rule concat_refseqs:
    input:
        results + "{subunit}_refseqs_gbff"
    output:
        results + "{subunit}_refseqs_proteomes.faa"
    shell:
        "cat {input}/ncbi_dataset/data/*/*.faa > {output}"

rule search_hzs_a_domain:
    input:
        proteomes=results + "hzs_bc_refseqs_proteomes.faa",
        hmm="../data/hzs-outside-anammox/PF18582.1.hmm"
    output:
        tbl=results + "hzs_bc_refseqs_proteomes.hzs_a_domanin.tsv",
        txt=results + "hzs_bc_refseqs_proteomes.hzs_a_domanin.txt"
    threads:
        12
    conda:
        "./envs/hmmer.yaml"
    shell:
        "hmmsearch --tblout {output.tbl} --noali -o {output.txt}  --cpu {threads} {input.hmm} {input.proteomes}"

rule extract_hzs_a_domain_proteins:
    input:
        results + "hzs_bc_refseqs_proteomes.hzs_a_domanin.tsv",
    output:
        results + "hzs_bc_refseqs_proteomes.hzs_a_domanin.filtered.tsv"
    shell:
        "python parse-hzs-a-hmmsearch.py {input} > {output}"

# Based on the HMM-search of the HZS-A middle domain, proteins with this domain
# was also identified in close proximity to the HZS-BC protein in Maribacter also.
# In the followin steps, those proteins are extracted and added to the HZS-A
# alignment and a new phylogeny is calculated.
#rule extract_maribacter_hzs_a:
#    input:
#
#    output:
#        results + "hzs_a.refseq.hq_anammox.maribacter.faa"

#    shell:

#rule add_maribacter_hzs_a:

#rule align_hzs_a:
#    input:
#        results + "hzs_a.refseq.hq_anammox.maribacter.faa"
#    output:
#        results + "hzs_a.refseq.hq_anammox.maribacter.aln"
#rule trim_hzs_a:
#
#rule iqtree_hzs_a:
# Plot Overview of HZS-BC-operon for anammox and sister group
#rule extract_hzs_bs
#    input:
#        proteins=
#    output:
#        results + "/"
#    shell:
#        "python  .py {input.genbanks} {input.proteins} > {output}"
# Extract all hits

rule extract_all_hits:
    input:
        results + "hzs.gtdb.w_taxa.tsv.gz"
    output:
        results + "{subunit}.gtdb.hits.faa"
    params:
        subunit = "{subunit}"
    shell:
        "python hzs-no-masking-extract-hits.py {input} {params.subunit} all all > {output}"

#rule run_interpro:
#    input:
#        results + "/{subunit}.gtdb.hits.faa"
#    output:
#        results + "/{subunit}.gtdb.hits.interpro.tsv"
#    shell:
#        "interproscan.sh -i {input} -o {output} -f TSV -cpu 12


# Create master table
rule merge_blast_and_structure:
    input:
        blast=results+"hzs.gtdb.w_taxa.tsv.gz",
        structure_alignment="../processed_data/hzs-structure-alignment-no-masking/alignment_summary.tsv",
        structure_map="../processed_data/hzs-structure-alignment-no-masking/{subunit}_alphafold_identifiers.tsv"
    output:
        results+"{subunit}.blast.structure.tsv"
    params:
        subunit="{subunit}"
    shell:
       "python merge_hzs_blast_structur_data.py {input.blast} {input.structure_map} {input.structure_alignment} {output} {params.subunit}"

rule add_MAG_info:
    input:
        table=results+"{subunit}.blast.structure.tsv",
        genome_info="../data/gtdb_representatives_database/gtdb_representatives.sample_gtdb.metadata.tsv"
    output:
        results + "{subunit}.blast.structure.genome.tsv"
    shell:
        "python add_MAG_info.py {input.table} {input.genome_info} {output}"

rule add_interpro_info:
    input:
        table=results + "{subunit}.blast.structure.genome.heme.tsv",
        interpro=results + "{subunit}.interproscan.tsv"
    output:
        results + "{subunit}.blast.structure.genome.heme.interpro.tsv"
    shell:
        "python add_interpro_info.py {input.interpro} {input.table} {output}"

rule add_heme_info:
    input:
        sequences=results+"{subunit}.gtdb.hits.faa",
        table=results+"{subunit}.blast.structure.genome.tsv",
    output:
        results + "{subunit}.blast.structure.genome.heme.tsv"
    conda:
        "./envs/biopython.yaml"
    shell:
        "python hzs-outside-anammox-cxxch.py {input.sequences} {input.table} {output}"

rule add_tree_annotation_column:
    input:
        local_taxonomy = "../data/local_taxonomy.tsv",
        local_seq = "../data/hzs-outside-anammox/{subunit}.faa",
        master=results + "{subunit}.blast.structure.genome.heme.tsv"
    output:
        master=results + "{subunit}.blast.structure.genome.heme.annotation.tsv"
    conda:
        "./envs/biopython.yaml"
    shell:
        "python merge_master_local_taxa.py {input.master} {input.local_taxonomy} {input.local_seq} {output}"
