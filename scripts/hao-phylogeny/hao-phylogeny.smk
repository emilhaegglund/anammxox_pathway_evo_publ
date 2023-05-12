"""
Workflow to reproduce the HAO-phylogeny.
"""
results = "../../processed-data/hao-phylogeny/"
data = "../../data/"
envs = "../envs/"

# Path to backup location
backup_dir = "/media/argos-emiha442/emiha442/1_projects/1_221123_anammox_pathway_evo/processed-data/hao-phylogeny/"

# Define genomes in the HQ-Dataset
HQ_GENOMES = ["GCA_000949635.1_ASM94963v1",
    "GCA_002443295.1_ASM244329v1",
    "GCA_007860005.1_ASM786000v1",
    "GCA_013112645.1_ASM1311264v1",
    "GCA_017347445.1_ASM1734744v1",
    "GCA_021650895.1_ASM2165089v1",
    "GCA_021650915.1_ASM2165091v1",
    "GCA_900232105.1_Kuenenia_stuttgartiensis_MBR1"]

rule all:
    input:
        # Blast search
        results + "hao-gtdb-hit-overlap.tsv",
        results + "hao-kust-identity.tsv",

        # Phylogenies
        results + "hao-gtdb-refseq-hq-anammox.remove-long-branch.trimal.fasttree",
        results + "hao-gtdb-refseq-hq-anammox.remove-long-branch.trimal.treefile",
        results + "hao-gtdb-refseq-hq-anammox.remove-long-branch.trimal-gt05.treefile",
        results + "hao-onr.trimal.treefile",
        results + "hao-onr-blast-annotation.tsv",
        results + "hao-onr-heme-clean-trimal.treefile",
        results + "hao-onr-heme-clean-trimal-gt05.treefile",
        results + "hao-hao_like.aln",

        # HAO6 cladde
        results + "hao6-clade.aln",
        results + "hao6-clade-all-vs-all-blast.tsv",
        #results + "hao-onr-cluster-trimal.treefile",

        # Uncomment following lines to backup
        #backup_dir + "hao-onr.trimal.treefile",
        backup_dir + "hao-onr-blast-annotation.tsv"

rule download_queries:
    "Download HAO-like proteins from Kuenenia"
    output:
        results + "query.faa"
    conda:
        envs + "entrez.yaml"
    shell:
        "efetch -id SOH04660.1,SOH06265.1,SOH05916.1,SOH05518.1,SOH05157.1,SOH04857.1,SOH03722.1,SOH04263.1,SOH05871.1,SOH05896.1 -db protein -format fasta > {output}"

rule blast_gtdb:
    "Search the GTDB-database for HAO-like proteins"
    input:
        results + "query.faa"
    output:
        results + "hao-gtdb.tsv.gz"
    params:
        db = data + "gtdb_representatives_database/gtdb_representatives.dmnd"
    conda:
        envs + "diamond.yaml"
    threads:
        30
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
    """Add the accession number to the HAO-hits, used to link it to the
       taxonomyin the next step."""
    input:
        blast = results + "hao-gtdb.tsv.gz",
        prot = data + "gtdb_representatives_database/prot2accession.tsv.gz"
    output:
        temp(results + "hao-gtdb-w-accessions.tsv.gz")
    shell:
        "python ../general/merge-blast-accessions.py {input.blast} {input.prot} {output}"

rule add_taxa:
    "Add taxonomic information to the HAO-hits"
    input:
        blast = results + "hao-gtdb-w-accessions.tsv.gz",
        gtdb = data + "gtdb_representatives_database/gtdb_representatives.sample_gtdb.metadata.tsv",
    output:
        results + "hao-gtdb-w-taxa.tsv.gz"
    shell:
        "python ../general/merge-blast-taxa.py {input.blast} {input.gtdb} {output}"

rule check_blast_overlap:
    "Calculate the Jaccard Index to check how similar "
    input:
        results + "hao-gtdb-w-taxa.tsv.gz"
    output:
        results + "hao-gtdb-hit-overlap.tsv",
    shell:
        "python check-diamond-hits-overlap.py {input} {output}"

rule hq_database:
    "Create a Diamond database for the high-quality proteomes"
    input:
        expand(data + "anammox-proteomes/{hq}_protein.faa", hq=HQ_GENOMES)
    output:
        results + "hq.dmnd"
    conda:
        envs + "diamond.yaml"
    shell:
        "cat {input} | diamond makedb --db {output} --in -"

rule blast_hq_dataset:
    "Search the anammox high-quality proteomes for HAO proteins"
    input:
        queries = results + "query.faa",
        db = results + "hq.dmnd"
    output:
        results + "hao-hq-anammox.tsv.gz"
    threads:
        12
    conda:
        envs + "diamond.yaml"
    shell:
        """
        diamond blastp \
            --db {input.db} \
            --query {input.queries} \
            --max-target-seqs 2000 \
            --outfmt 6 qseqid sseqid pident positive mismatch gapopen qlen slen length qstart qend sstart send evalue bitscore stitle full_sseq \
            --very-sensitive \
            --min-score 80 \
            --out {output} \
            --compress 1 \
            --threads {threads} \
        """

rule kust_blast_identity:
    """Extract the self hits of Kuenenia to get an estimate of
       how similar the HAO proteins inside Kuenenia are."""
    input:
        results + "hao-hq-anammox.tsv.gz"
    output:
        results + "hao-kust-identity.tsv"
    shell:
        "python hao-kust-blast.py {input} {output}"

rule extract_hq_hao_hits:
    "Extract HAO-hits from high-quality proteomes."
    input:
        results + "hao-hq-anammox.tsv.gz"
    output:
        results + "hao-hq-anammox.faa"
    shell:
        "python extract-hq-hao.py {input} > {output}"

rule extract_hao_gtdb_refseq:
    "Extract hits from search against GTDB"
    input:
        results + "hao-gtdb-w-taxa.tsv.gz"
    output:
        results + "hao-gtdb-refseq-no-anammox.faa"
    shell:
        """python ../general/extract-blast-hits.py --blast {input} \
            --fasta {output} \
            --qcover 0.5 \
            --scover 0.5 \
            --source refseq
        """

rule merge_anammox_hq_gtdb_refseq:
    input:
        anammox_hq = results + "hao-hq-anammox.faa",
        no_anammox = results + "hao-gtdb-refseq-no-anammox.faa"
    output:
        results + "hao-gtdb-refseq-hq-anammox.faa"
    shell:
        "cat {input.anammox_hq} {input.no_anammox} > {output}"

rule align:
    input:
        results + "hao-gtdb-refseq-hq-anammox.faa"
    output:
        results + "hao-gtdb-refseq-hq-anammox.aln"
    conda:
        envs + "mafft.yaml"
    threads:
        12
    shell:
        "mafft-linsi --thread {threads} {input} > {output}"

rule trim_alignment:
    input:
        results + "hao-gtdb-refseq-hq-anammox.aln"
    output:
        results + "hao-gtdb-refseq-hq-anammox.trimal.aln"
    conda:
        envs + "trimal.yaml"
    shell:
        "trimal -automated1 -in {input} -out {output}"

rule fasttree:
    input:
        results + "hao-gtdb-refseq-hq-anammox.trimal.aln"
    output:
        results + "hao-gtdb-refseq-hq-anammox.trimal.fasttree"
    conda:
        envs + "fasttree.yaml"
    shell:
        "fasttree -lg {input} > {output}"

rule remove_long_branch:
    input:
        sequences = results + "hao-gtdb-refseq-hq-anammox.faa",
        long_branches = data + "hao-phylogeny/long-branches.txt"
    output:
        results + "hao-gtdb-refseq-hq-anammox.remove-long-branch.faa",
    conda:
        envs + "biopython.yaml"
    shell:
        "python ../general/remove-taxa.py {input.sequences} {input.long_branches} {output}"

rule align_after_long_branch_remove_1:
    input:
        results + "hao-gtdb-refseq-hq-anammox.remove-long-branch.faa",
    output:
        results + "hao-gtdb-refseq-hq-anammox.remove-long-branch.aln",
    conda:
        envs + "mafft.yaml"
    threads:
        12
    shell:
        "mafft-linsi --thread {threads} {input} > {output}"

rule trim_after_long_branch_remove:
    input:
        results + "hao-gtdb-refseq-hq-anammox.remove-long-branch.aln",
    output:
        results + "hao-gtdb-refseq-hq-anammox.remove-long-branch.trimal.aln",
    conda:
        envs + "trimal.yaml"
    shell:
        "trimal -automated1 -in {input} -out {output}"

rule fasttree_after_long_branch_remove:
    input:
        results + "hao-gtdb-refseq-hq-anammox.remove-long-branch.trimal.aln",
    output:
        results + "hao-gtdb-refseq-hq-anammox.remove-long-branch.trimal.fasttree",
    conda:
        envs + "fasttree.yaml"
    shell:
        "fasttree -lg {input} > {output}"

rule iqtree:
    input:
        results + "hao-gtdb-refseq-hq-anammox.remove-long-branch.trimal.aln"
    output:
        results + "hao-gtdb-refseq-hq-anammox.remove-long-branch.trimal.treefile"
    params:
        pre = results + "hao-gtdb-refseq-hq-anammox.remove-long-branch.trimal"
    conda:
        envs + "iqtree.yaml"
    threads:
        12
    shell:
        "iqtree2 -s {input} -mset LG,WAG -B 1000 -alrt 1000 -nt {threads} -pre {params.pre} -redo"

rule gap_trim_after_long_branch_remove:
    input:
        results + "hao-gtdb-refseq-hq-anammox.remove-long-branch.aln",
    output:
        results + "hao-gtdb-refseq-hq-anammox.remove-long-branch.trimal-gt05.aln",
    conda:
        envs + "trimal.yaml"
    shell:
        "trimal -automated1 -in {input} -out {output}"

rule gap_trim_iqtree:
    input:
        results + "hao-gtdb-refseq-hq-anammox.remove-long-branch.trimal-gt05.aln"
    output:
        results + "hao-gtdb-refseq-hq-anammox.remove-long-branch.trimal-gt05.treefile"
    params:
        pre = results + "hao-gtdb-refseq-hq-anammox.remove-long-branch.trimal-gt05"
    conda:
        envs + "iqtree.yaml"
    threads:
        12
    shell:
        "iqtree2 -s {input} -mset LG,WAG -B 1000 -alrt 1000 -nt {threads} -pre {params.pre} -redo"
"""
Using MAD-rooting it has been suggested that HAO/HDH/IhOCC forms a separate
clade of Octaheme Cytochrome C proteins. Using this information we can root this
phylogeny with the Octaheme Nitrite reductase group.
"""
rule download_onr_outgroup_query:
    """
    Using the sequences as used in Soares et al. 2022 as queries to search
    the GTDB database.
    """
    output:
        results + "onr-outroup-query.faa"
    conda:
        envs + "entrez.yaml"
    shell:
        "efetch -id AGA31970.1,ADV18468.1 -db protein -format fasta > {output}"

rule blast_onr:
    input:
        results + "onr-outroup-query.faa"
    output:
        results + "onr-outgroup-gtdb.tsv.gz"
    params:
        db=data + "gtdb_representatives_database/gtdb_representatives.dmnd"
    conda:
        envs + "diamond.yaml"
    threads:
        30
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

rule onr_add_accessions:
    input:
        blast = results + "onr-outgroup-gtdb.tsv.gz",
        prot = data + "gtdb_representatives_database/prot2accession.tsv.gz"
    output:
        temp(results + "onr-outgroup-gtdb-w-accession.tsv.gz")
    shell:
        "python ../general/merge-blast-accessions.py {input.blast} {input.prot} {output}"

rule onr_add_taxa:
    input:
        blast = results + "onr-outgroup-gtdb-w-accession.tsv.gz",
        gtdb = data + "gtdb_representatives_database/gtdb_representatives.sample_gtdb.metadata.tsv",
    output:
        results + "onr-outgroup-gtdb-w-taxa.tsv.gz"
    shell:
        "python ../general/merge-blast-taxa.py {input.blast} {input.gtdb} {output}"

rule extract_onr_homologs:
    """
    Only use ADV18468.1, after only keep the best hit if detected by both
    queries only reports protein containing 5 hemes
    """
    input:
        results + "onr-outgroup-gtdb-w-taxa.tsv.gz"
    output:
        results + "onr-outgroup-gtdb-refseq.faa"
    shell:
        """
        python ../general/extract-blast-hits.py \
            --blast {input} \
            --queries ADV18468.1 \
            --scover 0.8 \
            --qcover 0.8 \
            --source refseq \
            --top 50 \
            --fasta {output}
        """

rule filter_cxxch:
    "Filter hits to contain eight hemes CxxC[H,K]"
    input:
        results + "onr-outgroup-gtdb-refseq.faa"
    output:
        table=results + "onr-outgroup-gtdb-refseq-heme-filtered.tsv",
        fasta=results + "onr-outgroup-gtdb-refseq-heme-filtered.faa"
    conda:
        envs + "biopython.yaml"
    shell:
        "python check-heme.py {input} {output.table} {output.fasta}"

rule merge_onr_hao:
    input:
        hao=results + "hao-gtdb-refseq-hq-anammox.remove-long-branch.faa",
        onr=results + "onr-outgroup-gtdb-refseq-heme-filtered.faa"
    output:
        results + "hao-onr.faa"
    shell:
        "cat {input.hao} {input.onr} > {output}"

rule align_onr_hao:
    input:
        results + "hao-onr.faa"
    output:
        results + "hao-onr.aln"
    conda:
        envs + "mafft.yaml"
    threads:
        12
    shell:
        "mafft-linsi --reorder --thread {threads} {input} > {output}"

rule trim_onr_hao:
    input:
        results + "hao-onr.aln"
    output:
        results + "hao-onr.trimal.aln"
    conda:
        envs + "trimal.yaml"
    shell:
        "trimal -automated1 -in {input} -out {output}"

rule fasttree_onr_hao:
    input:
        results + "hao-onr.trimal.aln"
    output:
        results + "hao-onr.trimal.fasttree"
    conda:
        envs + "fasttree.yaml"
    shell:
        "fasttree -lg {input} > {output}"

rule iqtree_onr_hao:
    input:
        results + "hao-onr.trimal.aln"
    output:
        results + "hao-onr.trimal.treefile"
    params:
        pre = results + "hao-onr.trimal"
    threads:
        12
    conda:
        envs + "iqtree.yaml"
    shell:
        "iqtree2 -s {input} -m MFP -mset LG,WAG -B 1000 -alrt 1000 -pre {params.pre} -nt {threads} -redo"

rule filter_hao_onr:
    "Filter one anammox clade with 11 heme groups and an early div group with 9 hemes"
    input:
        sequences = results + "hao-onr.faa",
        heme_clean = data + "hao-phylogeny/hao-onr-heme-filter.txt"
    output:
        results + "hao-onr-heme-clean.faa"
    conda:
        envs + "biopython.yaml"
    shell:
        "python ../general/remove-taxa.py {input.sequences} {input.heme_clean} {output}"

rule align_hao_onr_heme_clean:
    input:
        results + "hao-onr-heme-clean.faa"
    output:
        results + "hao-onr-heme-clean.aln"
    conda:
        envs + "mafft.yaml"
    threads:
        12
    shell:
        "mafft-linsi --thread {threads} --reorder {input} > {output}"

rule trim_hao_onr_heme_clean:
    input:
        results + "hao-onr-heme-clean.aln"
    output:
        results + "hao-onr-heme-clean-trimal.aln"
    conda:
        envs + "trimal.yaml"
    shell:
        "trimal -automated1 -in {input} -out {output}"

rule iqtree_hao_onr_heme_clean:
    input:
        results + "hao-onr-heme-clean-trimal.aln"
    output:
        results + "hao-onr-heme-clean-trimal.treefile"
    params:
        pre = results + "hao-onr-heme-clean-trimal"
    threads:
        12
    conda:
        envs + "iqtree.yaml"
    shell:
        "iqtree2 -s {input} -m MFP -mset LG,WAG -B 1000 -alrt 1000 -pre {params.pre} -nt {threads} -redo"

rule gap_trim_hao_onr_heme_clean:
    input:
        results + "hao-onr-heme-clean.aln"
    output:
        results + "hao-onr-heme-clean-trimal-gt05.aln"
    conda:
        envs + "trimal.yaml"
    shell:
        "trimal -gt 0.5 -in {input} -out {output}"

rule iqtree_gap_trim_hao_onr_heme_clean:
    input:
        results + "hao-onr-heme-clean-trimal-gt05.aln"
    output:
        results + "hao-onr-heme-clean-trimal-gt05.treefile"
    params:
        pre = results + "hao-onr-heme-clean-trimal-gt05"
    threads:
        12
    conda:
        envs + "iqtree.yaml"
    shell:
        "iqtree2 -s {input} -m MFP -mset LG,WAG -B 1000 -alrt 1000 -pre {params.pre} -nt {threads} -redo"

rule cluster_hao_onr:
    """
    In this step the proteins are clustered at 80% identity. This reduces the
    number of sequences from 815 to 418 sequences.
    """
    input:
        results + "hao-onr.faa"
    output:
        results + "hao-onr-cluster.faa"
    conda:
        envs + "cd-hit.yaml"
    shell:
        "cd-hit -i {input} -o {output} -c 0.8"

rule align_cluster_hao_onr:
    """
    This alignment contains 2602 sites.
    """
    input:
        results + "hao-onr-cluster.faa"
    output:
        results + "hao-onr-cluster.aln"
    conda:
        envs + "mafft.yaml"
    threads:
        12
    shell:
        "mafft-linsi --reorder --thread {threads} {input} > {output}"

rule trim_clustered_hao_onr:
    "After trimming there are 254 sites in the alignment"
    input:
        results + "hao-onr-cluster.aln"
    output:
        results + "hao-onr-cluster-trimal.aln"
    conda:
        envs + "trimal.yaml"
    shell:
        "trimal -automated1 -in {input} -out {output}"

rule iqtree_clustered_hao_onr:
    input:
        results + "hao-onr-cluster-trimal.aln"
    output:
        results + "hao-onr-cluster-trimal.treefile"
    params:
        pre = results + "hao-onr-cluster-trimal"
    threads:
        12
    conda:
        envs + "iqtree.yaml"
    shell:
        "iqtree2 -s {input} -m MFP -mset LG,WAG -B 1000 -alrt 1000 -pre {params.pre} -nt {threads} -redo"

rule blast_table:
    input:
        hao = results + "hao-gtdb-w-taxa.tsv.gz",
        onr = results + "onr-outgroup-gtdb-w-taxa.tsv.gz"
    output:
        results + "hao-onr-blast.tsv"
    shell:
        "python create-hao-blast-table.py {input} {output}"

rule add_tree_annotation_column:
    input:
        local_taxonomy = data + "local-taxonomy.tsv",
        local_seq = results + "hao-hq-anammox.faa",
        master = results + "hao-onr-blast.tsv"
    output:
        results + "hao-onr-blast-annotation.tsv"
    conda:
        envs + "biopython.yaml"
    shell:
        "python ../general/merge-blast-table-local-taxa.py {input.master} {input.local_taxonomy} {input.local_seq} {output}"

rule extract_hao_hao_like_sequences:
    input:
        accessions = data + "hao-phylogeny/hao-hao_like-accessions.txt",
        sequences = results + "hao-onr-heme-clean.faa"
    output:
        results + "hao-hao_like.faa"
    conda:
        envs + "seqtk.yaml"
    shell:
        "seqtk subseq {input.sequences} {input.accessions} > {output}"

rule aling_hao_hao_like:
    input:
        results + "hao-hao_like.faa"
    output:
        results + "hao-hao_like.aln"
    conda:
        envs + "mafft.yaml"
    threads:
        12
    shell:
        "mafft-linsi --thread {threads} {input} > {output}"

rule extract_hao6_group:
    input:
        accessions = data + "hao-phylogeny/hao-clades.tsv",
        sequences = results + "hao-onr-heme-clean.faa"
    output:
        results + "hao6-clade.faa"
    conda:
        envs + "seqtk.yaml"
    shell:
        "grep 'HAO6' {input.accessions} | awk -F'\\t' '{{ print $1  }}' | seqtk subseq {input.sequences} - > {output}"

rule align_hao6_clade:
    input:
        results + "hao6-clade.faa"
    output:
        results + "hao6-clade.aln"
    conda:
        envs + "mafft.yaml"
    shell:
        "mafft-linsi {input} > {output}"

rule hao6_blast:
    input:
        db = results + "hq.dmnd",
        queries = results + "hao6-clade.faa"
    output:
        results + "hao6-clade-all-vs-all-blast.tsv"
    conda:
        envs + "diamond.yaml"
    threads:
        12
    shell:
        """
        diamond blastp \
            --db {input.db} \
            --query {input.queries} \
            --max-target-seqs 2000 \
            --outfmt 6 qseqid sseqid pident positive mismatch gapopen qlen slen length qstart qend sstart send evalue bitscore stitle full_sseq \
            --very-sensitive \
            --min-score 80 \
            --out {output} \
            --threads {threads} \
        """

rule backup:
    input:
        hao_onr_tree = results + "hao-onr.trimal.treefile",
        tree_annotation = results + "hao-onr-blast-annotation.tsv"
    output:
        backup_dir + "hao-onr-blast-annotation.tsv",
        backup_dir + "hao-onr.trimal.treefile"
    params:
        backup_dir = backup_dir
    shell:
        "cp {input.hao_onr_tree} {input.tree_annotation} {params.backup_dir}"
