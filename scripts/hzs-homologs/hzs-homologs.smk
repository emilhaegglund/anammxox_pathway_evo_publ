backup_dir = "/media/argos-emiha442/emiha442/1_projects/1_221123_anammox_pathway_evo/processed_data/hzs-homologs/"
results = "../processed-data/hzs-homologs/"
data = "../../data/"
envs = "../envs/"
bin_dir = "../../bin"  # Path to directory containing software not in Bioconda
SUBUNITS = ["hzs_a", "hzs_bc"]
rule all:
    input:
        # Summary of blast search
        results + "hzs-blast-summary.svg",

        # Summary of structure alignment
        os.path.join(backup_dir, "structure-alignment/alignment-summary.pdf"),
        os.path.join(backup_dir, "structure-alignment/alignment-summary.tsv")
        
        # Phylogenies
        expand(results + "{subunit}-gtdb-refseq-long-branch-remove.trimal.rooted.nwk", subunit=SUBUNITS),
        #expand(results + "{subunit}.blast.structure.genome.heme.annotation.tsv", subunit=SUBUNITS),
        #expand(results + "{subunit}.gtdb.hits.faa", subunit=SUBUNITS),
        #results + "hzs_bc.sister_group.anammox.aln",
        #results + "hzs_bc_refseqs_proteomes.faa",
        #results + "hzs_bc_refseqs_proteomes.hzs_a_domanin.filtered.tsv",
        #results + "hzs_a.refseq.hq_anammox.maribacter.faa"

rule download_queries:
    "Download HZS-B, HZS-C, and HZS-A from Kuenenia stuttgartiensis and HZS-BC and HZS-A from Scalinduae japonica"
    output:
        results + "query.faa"
    conda:
        envs + "entrez.yaml"
    shell:
        "efetch -id SOH05198.1,SOH05199.1,SOH05200.1,GAX62881.1,GAX62882.1 -db protein -format fasta > {output}"

rule blast_gtdb:
    "Search for homologous proteins among the representative species in GTDB v207"
    input:
        results + "query.faa"
    output:
        results + "hzs-gtdb.tsv.gz"
    params:
        db = data + "gtdb_representatives_database/gtdb_representatives.dmnd"
    conda:
        envs + "diamond.yaml"
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
        blast = results + "hzs-gtdb.tsv.gz",
        prot = data + "gtdb_representatives_database/prot2accession.tsv.gz"
    output:
        temp(results + "hzs-gtdb-w-accessions.tsv.gz")
    shell:
        "python ../general/merge-blast-accessions.py {input.blast} {input.prot} {output}"

rule add_taxa:
    "Add taxonomic information for each of the identified proteins"
    input:
        blast = results + "hzs-gtdb-w-accessions.tsv.gz",
        gtdb = data + "gtdb_representatives_database/gtdb_representatives.sample_gtdb.metadata.tsv",
    output:
        results + "hzs-gtdb-w-taxa.tsv.gz"
    shell:
        "python ../general/merge-blast-taxa.py {input.blast} {input.gtdb} {output}"

rule plot_summary_blast_results:
    "Plot a summary of the blast results. Presented in supplementary figure X"
    input:
        results + "hzs-gtdb-w-taxa.tsv.gz"
    output:
        results + "hzs-blast-summary.svg"
    shell:
        "python plot-blast-results-distribution.py {input} {output}"

# Run pairwise structure alignment between query and hits
rule find_alpha_fold_structures:
    """
    Parse uniprot and collect links to AlphaFold structures for the hits
    which have structures available.
    """
    input:
        blast_df = results + "hzs.gtdb.w_taxa.tsv.gz"
    output:
        results + "structure-alignment/{subunit}-alphafold-identifiers.tsv"
    params:
        subunit = "{subunit}"
    shell:
        "python search-structures-uniprot.py {input} {output} {params.subunit}"

rule download_hzs_a_structure:
    """
    Download Scalinduae 
    """
    output:
        results + "structure-alignment/hzs_a-af-structure.pdb"
    shell:
        "wget -O {output} https://alphafold.ebi.ac.uk/files/AF-A0A286U438-F1-model_v4.pdb"

rule download_hzs_bc_structure:
    """
    Download Scalinduae 
    """
    output:
        results + "structure-alignment/hzs_bc-af-structure.pdb"
    shell:
        "wget -O {output} https://alphafold.ebi.ac.uk/files/AF-A0A286U486-F1-model_v4.pdb"

checkpoint download_homolog_hzs_bc_structures:
    input:
        results + "structure-alignment/hzs_bc-alphafold-identifiers.tsv"
    output:
        directory(results + "structure-alignment/hzs_bc-structures")
    shell:
        """
        awk -F'\\t' '{{ print $2 }}' {input} | grep -v 'alphafold_url' | sed '/^$/d' | wget -i - -P {output}
        """

def get_hzs_bc_structures(wildcards):
    ck_output = checkpoints.download_homolog_hzs_bc_structures.get(**wildcards).output[0]
    homologs, = glob_wildcards(os.path.join(ck_output, "{homologs}.pdb"))
    return os.path.join(results, "structure-alignment", "hzs_bc-structures",
                               "{homologs}.pdb")

rule structure_alignment_hzs_bc:
    input:
        ref = results + "structure-alignment/hzs_bc-af-structure.pdb",
        target = get_hzs_bc_structures
    output:
        results + "structure-alignment/hzs_bc-structure-alignment/{homologs}.txt"
    conda:
        envs + "tmalign.yaml"
    shell:
        """
        TMalign {input.ref} {input.target} -a > {output}
        """

checkpoint download_homolog_hzs_a_structures:
    input:
        results + "structure-alignment/hzs_a-alphafold-identifiers.tsv"
    output:
        directory(results + "structure-alignment/hzs_a-structures")
    shell:
        """
        awk -F'\\t' '{{ print $2 }}' {input} | grep -v 'alphafold_url' | sed '/^$/d' | wget -i - -P {output}
        """

def get_hzs_a_structures(wildcards):
    ck_output = checkpoints.download_homolog_hzs_a_structures.get(**wildcards).output[0]
    homologs, = glob_wildcards(os.path.join(ck_output, "{homologs}.pdb"))
    return os.path.join(results, "structure-alignment", "hzs_a-structures",
                               "{homologs}.pdb")

rule structure_alignment_hzs_a:
    input:
        ref = results + "structure-alignment/hzs_a-af-structure.pdb",
        target= get_hzs_a_structures
    output:
        results + "structure-alignment/hzs_a-structure-alignment/{homologs}.txt"
    conda:
        envs + "tmalign.yaml"
    shell:
        """
        TMalign {input.ref} {input.target} -a > {output}
        """

def get_hzs_bc_structure_alignment(wildcards):
    ck_output = checkpoints.download_homolog_hzs_bc_structures.get(**wildcards).output[0]
    homologs, = glob_wildcards(os.path.join(ck_output, "{homologs}.pdb"))
    return expand(os.path.join(results, "structure-alignment", "hzs_bc-structure-alignment", "{homologs}.txt"), homologs=homologs)

def get_hzs_a_structure_alignment(wildcards):
    ck_output = checkpoints.download_homolog_hzs_a_structures.get(**wildcards).output[0]
    homologs, = glob_wildcards(os.path.join(ck_output, "{homologs}.pdb"))
    return expand(os.path.join(results, "structure-alignment", "hzs_a-structure-alignment", "{homologs}.txt"), homologs=homologs)

rule summarize_results:
    input:
        get_hzs_bc_structure_alignment,
        get_hzs_a_structure_alignment
    output:
        results + "structure-alignment/alignment-summary.tsv"
    shell:
        "python parse-structure-alignments.py {input} {output}"

rule plot_results:
    input:
        results + "structure-alignment/alignment-summary.tsv"
    output:
        results + "structure-alignment/alignment-summary.pdf"
    shell:
        "python plot-structure-alignment-summary.py {input} {output}"

rule backup_results:
    input:
        results + "structure-alignment/alignment_summary.tsv"
    output:
        backup_dir +  "structure-alignemt/alignment-summary.tsv"
    shell:
        "cp {input} {output}"

rule backup_figure:
    input:
        results + "structure-alignment/alignment-summary.pdf"
    output:
        backup_dir + "structure-alignment/alignment-summary.pdf")
    shell:
        "cp {input} {output}"

# Prepare phylogeny from RefSeq homologs
rule extract_hzs_a_refseq_sequences:
    """Extract protein sequence for all hits to HZS-A for entries derived from
    RefSeq"""
    input:
        blast = results + "hzs-gtdb-w-taxa.tsv.gz"
    output:
        results + "hzs_a-gtdb-refseq-no-anammox.faa"
    shell:
        """python ../general/extract-blast-hits.py \
            --blast {input.blast} \
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
        blast = results + "hzs-gtdb-w-taxa.tsv.gz"
    output:
        results + "hzs_bc-gtdb-refseq-no-anammox.faa"
    shell:
        """python extract-blast-hits.py \
            --blast {input.blast} \
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
        no_anammox = results + "{subunit}-gtdb-refseq-no-anammox.faa",
        anammox = data + "hzs-homologs/{subunit}.faa"
    output:
        results + "{subunit}-gtdb-refseq.faa"
    shell:
        "cat {input.anammox} {input.no_anammox} > {output}"

rule align_refseqs:
    "Align the sequences using mafft."
    input:
        results + "{subunit}-gtdb-refseq.faa"
    output:
        results + "{subunit}-gtdb-refseq.aln"
    conda:
        envs + "mafft.yaml"
    threads:
        12
    shell:
        "mafft-linsi --reorder --thread {threads} {input} > {output}"

rule trim_refseqs:
    "Trim the aligment using TrimAl."
    input:
        results + "{subunit}-gtdb-refseq.aln"
    output:
        results + "{subunit}-gtdb-refseq.trimal.aln"
    conda:
        envs + "trimal.yaml"
    threads:
        12
    shell:
        "trimal -automated1 -keepheader -in {input} -out {output}"

rule fasttree_refseqs:
    "Run FastTree to get first overview of data."
    input:
        results + "{subunit}-gtdb-refseq.trimal.aln"
    output:
        results + "{subunit}-gtdb-refseq.trimal.fasttree"
    conda:
        envs + "fasttree.yaml"
    shell:
        "fasttree -lg {input} > {output}"

rule remove_gappy_sequences_refseqs:
    input:
        aln = results + "{subunit}-gtdb-refseq.trimal.aln",
        seq = results + "{subunit}-gtdb-refseq.faa"
    output:
        results + "{subunit}-gtdb-refseq-gap-remove.faa"
    params:
        gap_frac = 0.5
    conda:
        envs + "biopython.yaml"
    shell:
        "python ../general/remove-gappy-sequences.py {input.aln} {input.seq} {params.gap_frac} {output}"

rule align_refseqs_after_gap_remove:
    input:
        results + "{subunit}-gtdb-refseq-gap-remove.faa"
    output:
        results + "{subunit}-gtdb-refseq-gap-remove.aln"
    conda:
        envs + "mafft.yaml"
    threads:
        12
    shell:
        "mafft-linsi --reorder --thread {threads} {input} > {output}"

rule trim_refseqs_after_gap_remove:
    input:
        results + "{subunit}-gtdb-refseq-gap-remove.aln"
    output:
        results + "{subunit}-gtdb-refseq-gap-remove.trimal.aln"
    conda:
        envs + "trimal.yaml"
    threads:
        12
    shell:
        "trimal -automated1 -keepheader -in {input} -out {output}"

rule fasttree_refseqs_after_gap_remove:
    input:
        results + "{subunit}-gtdb-refseq-gap-remove.trimal.aln"
    output:
        results + "{subunit}-gtdb-refseq-gap-remove.trimal.fasttree"
    conda:
        envs + "fasttree.yaml"
    shell:
        "fasttree -lg {input} > {output}"

rule remove_long_branches:
    input:
        sequences = results + "{subunit}-gtdb-refseq-gap-remove.faa",
        accessions = data + "hzs-homologs/{subunit}-long-branches.txt"
    output:
        results + "{subunit}-gtdb-refseq-long-branch-remove.faa"
    conda:
        envs + "biopython.yaml"
    shell:
        "python ../general/remove-taxa.py {input.sequences} {input.accessions} {output}"

rule align_refseqs_after_remove_long_branches:
    input:
        results + "{subunit}-gtdb-refseq-long-branch-remove.faa"
    output:
        results + "{subunit}-gtdb-refseq-long-branch-remove.aln"
    conda:
        envs + "mafft.yaml"
    threads:
        12
    shell:
        "mafft-linsi --reorder --thread {threads} {input} > {output}"

rule trim_refseqs_after_remove_long_branches:
    input:
        results + "{subunit}-gtdb-refseq-long-branch-remove.aln"
    output:
        results + "{subunit}-gtdb-refseq-long-branch-remove.trimal.aln"
    conda:
        envs + "trimal.yaml"
    threads:
        12
    shell:
        "trimal -automated1 -keepheader -in {input} -out {output}"

rule fasttree_refseqs_after_remove_long_branches:
    input:
        results + "{subunit}-gtdb-refseq-long-branch-remove.trimal.aln"
    output:
        results + "{subunit}-gtdb-refseq-long-branch-remove.trimal.fasttree"
    conda:
        envs + "fasttree.yaml"
    shell:
        "fasttree -lg {input} > {output}"

rule iqtree_refseqs:
    input:
        results + "{subunit}-gtdb-refseq-long-branch-remove.trimal.aln"
    output:
        results + "{subunit}-gtdb-refseq-long-branch-remove.trimal.treefile"
    conda:
        envs + "iqtree.yaml"
    params:
        pre = results + "{subunit}-gtdb-refseq-long-branch-remove.trimal"
    threads:
        12
    shell:
        "iqtree2 -s {input} -mset LG,WAG -alrt 1000 -B 1000 -nt {threads} -pre {params.pre} -redo"

rule mad_root:
    input:
        results + "{subunit}-gtdb-refseq-long-branch-remove.trimal.treefile"
    output:
        results + "{subunit}-gtdb-refseq-long-branch-remove.trimal.rooted.nwk"
    shell:
        bin_dir + "MADroot/madRoot {input} > {output}"

# Align only Brocadiae and sister group
rule extract_sister_group:
    input:
        accessions = data + "hzs-homologs/hzs_bc-sister-group.txt",
        sequences = results + "hzs_bc-gtdb-refseq.faa"
    output:
        results + "hzs_bc-sister-group.faa"
    conda:
        envs + "seqtk.yaml"
    shell:
        "seqtk subseq {input.sequences} {input.accessions} > {output}"

rule merge_sister_group_anammox:
    input:
        sister_group = results + "hzs_bc-sister-group.faa",
        anammox = data + "hzs-homologs/hzs_bc.faa"
    output:
        results + "hzs_bc-sister-group-anammox.faa"
    shell:
        "cat {input.sister_group} {input.anammox} > {output}"

rule align_anammox_sister_group:
    input:
        results + "hzs_bc-sister-group-anammox.faa"
    output:
        results + "hzs_bc-sister-group-anammox.aln"
    conda:
        envs + "mafft.yaml"
    shell:
        "mafft-linsi --reorder {input} > {output}"


# Download GBFF for comparative analysis of the Refseqs. The following steps
# download proteomes and GBK files for all proteomes in RefSeq with a hit to the
# HZS-BC subunit. Since HZS-A is less conserved compared to HZS-BC we search the
# proteomes for the HZS-A middle domain using hmmsearch. The GBKs is also used
# to plot the genome structure around the HZS-BC homologs.
rule extract_refseq_accessions:
    input:
        results + "{subunit}-blast-structure-genome.tsv"
    output:
        results + "{subunit}-gtdb-refseq-accessions.txt"
    params:
        subunit = "{subunit}"
    shell:
        "python extract-refseq-accession.py {input} {params.subunit} > {output}"

# The following was run outside the Snakemake workflow.
#rule download_refseq_gbff:
#    input:
#        results + "{subunit}-gtdb-refseq-accessions.txt"
#    output:
#        results + "{subunit}-gtdb-refseq.zip"
#    conda:
#        envs + "ncbi-datasets.yaml"
#    shell:
#        "datasets download genome accession --inputfile {input} --filename {output} --include gbff,protein"

rule unzip_refseq_gbff:
    input:
        results + "{subunit}-gtdb-refseq.zip"
    output:
        directory(results + "{subunit}-gtdb-refseq-genome-data")
    shell:
        "unzip {input} -d {output}"

rule concat_refseqs:
    input:
        results + "{subunit}-gtdb-refseq-genome-data"
    output:
        results + "{subunit}-gtdb-refseq-proteomes.faa"
    shell:
        "cat {input}/ncbi_dataset/data/*/*.faa > {output}"

rule search_hzs_a_domain:
    input:
        proteomes = results + "hzs_bc-gtdb-refseq-proteomes.faa",
        hmm = data + "hzs-homologs/PF18582.1.hmm"
    output:
        tbl = results + "hzs_bc-gtdb-refseq-proteomes-hzs_a-domanin.tsv",
        txt = results + "hzs_bc-gtdb-refseq_proteomes-hzs_a-domanin.txt"
    threads:
        12
    conda:
        envs + "hmmer.yaml"
    shell:
        "hmmsearch --tblout {output.tbl} --noali -o {output.txt}  --cpu {threads} {input.hmm} {input.proteomes}"

rule extract_hzs_a_domain_proteins:
    input:
        results + "hzs_bc-gtdb-refseq-proteomes-hzs_a-domanin.tsv",
    output:
        results + "hzs_bc-gtdb-refseq_proteomes-hzs_a-domanin-filtered.tsv"
    shell:
        "python parse-hzs-a-hmmsearch.py {input} > {output}"

# Based on the HMM-search of the HZS-A middle domain, proteins with this domain
# was also identified in close proximity to the HZS-BC protein in Maribacter also.
# In the followin steps, those proteins are extracted and added to the HZS-A
# alignment and a new phylogeny is calculated.
rule extract_maribacter_hzs_a:
    output:
        results + "hzs_a-maribacter.faa"
    conda:
        envs + "entrez.yaml"
    shell:
        "efetch -id WP_188243417.1,WP_188314372.1,WP_109649145.1,WP_154366457.1,WP_163345621.1 -db protein -format fasta > {output}"

rule add_maribacter_hzs_a:
    input:
        maribacter = results+"hzs_a-maribacter.faa",
        anammox = results + "hzs_a-gtdb-refseq-long-branch-remove.faa"
    output:
        results + "hzs_a-gtdb-refseq-maribacter.faa"
    shell:
        "cat {input.maribacter} {input.anammox} > {output}"

rule align_maribacter_hzs_a:
    input:
        results + "hzs_a-gtdb-refseq-maribacter.faa"
    output:
        results + "hzs_a-gtdb-refseq-maribacter.aln"
    conda:
        envs + "mafft.yaml"
    threads:
        12
    shell:
        "mafft-linsi --reorder --thread {threads} {input} > {output}"

rule trim_maribacter_hzs_a:
    input:
        results + "hzs_a-gtdb-refseq-maribacter.aln"
    output:
        results + "hzs_a-gtdb-refseq-maribacter-trimal.aln"
    conda:
        envs + "trimal.yaml"
    shell:
        "trimal -automated1 -in {input} -out {output}"

rule iqtree_maribacterr_hzs_a:
    input:
        results + "hzs_a-gtdb-refseq-maribacter-trimal.aln"
    output:
        results + "hzs_a-gtdb-refseq-maribacter-trimal.treefile"
    conda:
        envs + "iqtree.yaml"
    params:
        pre = results + "{subunit}-gtdb-refseq-maribacter-trimal"
    threads:
        12
    shell:
        "iqtree2 -s {input} -mset LG,WAG -alrt 1000 -B 1000 -nt {threads} -pre {params.pre} -redo"

rule mad_root_maribacter:
    input:
        results + "hzs_a-gtdb-refseq-maribacter-trimal.treefile"
    output:
        results + "hzs_a-gtdb-refseq-maribacter-trimal-rooted.nwk"
    shell:
        bin_dir + "MADroot/madRoot {input} > {output}"

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
        blast = results + "hzs-gtdb-w-taxa.tsv.gz",
        structure_alignment = results + "structure-alignment/alignment-summary.tsv",
        structure_map = results + "structure-alignment/{subunit}-alphafold-identifiers.tsv"
    output:
        results + "{subunit}-blast-structure.tsv"
    params:
        subunit = "{subunit}"
    shell:
       "python merge-hzs-blast-structure-alignment.py {input.blast} {input.structure_map} {input.structure_alignment} {output} {params.subunit}"

rule add_MAG_info:
    input:
        table = results + "{subunit}-blast-structure.tsv",
        genome_info = data + "gtdb_representatives_database/gtdb_representatives.sample_gtdb.metadata.tsv"
    output:
        results + "{subunit}-blast-structure-genome.tsv"
    shell:
        "python add-genome-metadata.py {input.table} {input.genome_info} {output}"

rule add_heme_info:
    input:
        sequences = results + "{subunit}-gtdb-refseq.faa",
        table = results + "{subunit}-blast-structure-genome.tsv",
    output:
        results + "{subunit}-blast-structure-genome-heme.tsv"
    conda:
        envs + "biopython.yaml"
    shell:
        "python count-cxxch.py {input.sequences} {input.table} {output}"

rule add_interpro_info:
    input:
        table = results + "{subunit}-blast-structure-genome-heme.tsv",
        interpro = results + "{subunit}-interproscan.tsv"
    output:
        results + "{subunit}-blast-structure-genome-heme-interpro.tsv"
    shell:
        "python add-interpro-info.py {input.interpro} {input.table} {output}"

rule add_tree_annotation_column:
    input:
        local_taxonomy = data + "local_taxonomy.tsv",
        local_seq = data + "hzs-homologs/{subunit}.faa",
        master = results + "{subunit}-blast-structure-genome-heme.tsv"
    output:
        results + "{subunit}-blast-structure-genome-heme-annotation.tsv"
    conda:
        envs + "biopython.yaml"
    shell:
        "python merge-master-table-and-local-taxa.py {input.master} {input.local_taxonomy} {input.local_seq} {output}"