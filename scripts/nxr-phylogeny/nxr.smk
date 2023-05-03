backup_dir = "/media/argos-emiha442/emiha442/1_projects/1_221123_anammox_pathway_evo/processed_data/nxr-phylogeny/"
results = "../../processed-data/nxr-phylogeny/"
data = "../../data/"
envs = "../envs/"

subunits = [
    "nxr_a",
    "nxr_b",
    "nxr_c",
    "nxr_d",
    "nxr_hypo",
    "nxr_uspa",
]
rule all:
    input:
        expand(results + "{subunit}-gtdb-refseq-hq-anammox-clustered-trimal.fasttree", subunit=subunits),
        expand(results + "{subunit}-gtdb-refseq-hq-anammox-clustered-trimal.treefile", subunit=subunits),
        expand(results + "{subunit}-gtdb-refseq-clade-outgroup-trimal.treefile", subunit=["nxr_a"]),
        expand(results + "{subunit}-phylo-annotation.tsv", subunit=subunits),

        # MAGs
        #expand(results + "{subunit}-gtdb-clustered-trimal.fasttree", subunit=["nxr_c", "nxr_d", "nxr_hypo"]),
        #expand(results + "{subunit}-gtdb.tsv", subunit=["nxr_c", "nxr_d", "nxr_hypo"])



rule download_queries:
    output:
        results + "nxr-query.faa"
    conda:
        envs + "entrez.yaml"
    shell:
        "efetch -id SOH03954.1,SOH03955.1,SOH03956.1,SOH03957.1,SOH03958.1,SOH03960.1 -db protein -format fasta > {output}"

rule blast_gtdb:
    input:
        results + "nxr-query.faa"
    output:
        results + "nxr-gtdb.tsv.gz"
    params:
        db="../../data/gtdb_representatives_database/gtdb_representatives.dmnd"
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
    input:
        blast = results + "nxr-gtdb.tsv.gz",
        prot = data + "gtdb_representatives_database/prot2accession.tsv.gz"
    output:
        temp(results + "nxr-gtdb-w-accessions.tsv.gz")
    shell:
        "python ../general/merge-blast-accessions.py {input.blast} {input.prot} {output}"

rule add_taxa:
    input:
        blast = results + "nxr-gtdb-w-accessions.tsv.gz",
        gtdb = data + "gtdb_representatives_database/gtdb_representatives.sample_gtdb.metadata.tsv",
    output:
        results + "nxr-gtdb-w-taxa.tsv.gz"
    shell:
        "python ../general/merge-blast-taxa.py {input.blast} {input.gtdb} {output}"

rule extract_nxr_a_gtdb_refseq:
    input:
        results + "nxr-gtdb-w-taxa.tsv.gz"
    output:
        sequences=results + "nxr_a-gtdb-refseq-no-anammox.faa",
        table=results + "nxr_a-gtdb-refseq-no-anammox.tsv"
    shell:
        """
        python ../general/extract-blast-hits.py \
            --blast {input} \
            --fasta {output.sequences} \
            --table {output.table} \
            --queries SOH03954.1 \
            --qcover 0.5 \
            --scover 0.5 \
            --source refseq
        """

rule extract_nxr_b_gtdb_refseq:
    input:
        results + "nxr-gtdb-w-taxa.tsv.gz",
    output:
        sequences=results + "nxr_b-gtdb-refseq-no-anammox.faa",
        table=results + "nxr_b-gtdb-refseq-no-anammox.tsv"
    shell:
        """
        python ../general/extract-blast-hits.py \
            --blast {input} \
            --fasta {output.sequences} \
            --table {output.table} \
            --queries SOH03957.1 \
            --qcover 0.5 \
            --scover 0.5 \
            --source refseq
        """

rule extract_nxr_c_gtdb_refseq:
    input:
        results + "nxr-gtdb-w-taxa.tsv.gz"
    output:
        sequences=results + "nxr_c-gtdb-refseq-no-anammox.faa",
        table=results + "nxr_c-gtdb-refseq-no-anammox.tsv"
    shell:
        """
        python ../general/extract-blast-hits.py \
            --blast {input} \
            --fasta {output.sequences} \
            --table {output.table} \
            --queries SOH03958.1 \
            --qcover 0.5 \
            --scover 0.5 \
            --source refseq
        """

rule extract_nxr_d_gtdb_refseq:
    input:
        results + "nxr-gtdb-w-taxa.tsv.gz"
    output:
        sequences=results + "nxr_d-gtdb-refseq-no-anammox.faa",
        table=results + "nxr_d-gtdb-refseq-no-anammox.tsv"
    shell:
        """
        python ../general/extract-blast-hits.py \
            --blast {input} \
            --fasta {output.sequences} \
            --table {output.table} \
            --queries SOH03956.1 \
            --qcover 0.5 \
            --scover 0.5 \
            --source refseq
        """

rule extract_nxr_hypo_gtdb_refseq:
    input:
        results + "nxr-gtdb-w-taxa.tsv.gz"
    output:
        sequences=results + "nxr_hypo-gtdb-refseq-no-anammox.faa",
        table=results + "nxr_hypo-gtdb-refseq-no-anammox.tsv"
    shell:
        """
        python ../general/extract-blast-hits.py \
            --blast {input} \
            --fasta {output.sequences} \
            --table {output.table} \
            --queries SOH03955.1 \
            --qcover 0.5 \
            --scover 0.5 \
            --source refseq
        """

rule extract_nxr_uspa_gtdb_refseq:
    input:
        results + "nxr-gtdb-w-taxa.tsv.gz"
    output:
        sequences=results + "nxr_uspa-gtdb-refseq-no-anammox.faa",
        table=results + "nxr_uspa-gtdb-refseq-no-anammox.tsv"
    shell:
        """
        python ../general/extract-blast-hits.py \
            --blast {input} \
            --fasta {output.sequences} \
            --table {output.table} \
            --queries SOH03960.1 \
            --qcover 0.5 \
            --scover 0.5 \
            --source refseq
        """

rule merge_anammox_no_anammox:
    input:
        no_anammox = results + "{subunit}-gtdb-refseq-no-anammox.faa",
        anammox = data + "nxr-phylogeny/{subunit}-hq-anammox.faa"
    output:
        results + "{subunit}-gtdb-refseq-hq-anammox.faa"
    shell:
        "cat {input.no_anammox} {input.anammox} > {output}"

rule cluster_nxr_subunits:
    input:
        results + "{subunit}-gtdb-refseq-hq-anammox.faa"
    output:
        sequences = results + "{subunit}-gtdb-refseq-hq-anammox-clustered.faa",
        clusters = results + "{subunit}-gtdb-refseq-hq-anammox-clustered.faa.clstr"
    conda:
        "../envs/cd-hit.yaml"
    shell:
        "cd-hit -i {input} -o {output.sequences} -c 0.9"

rule align:
    input:
        results + "{subunit}-gtdb-refseq-hq-anammox-clustered.faa"
    output:
        results + "{subunit}-gtdb-refseq-hq-anammox-clustered.aln"
    conda:
        "../envs/mafft.yaml"
    threads:
        12
    shell:
        "mafft-linsi --thread {threads} --reorder {input} > {output}"

rule trim_alignment:
    input:
        results + "{subunit}-gtdb-refseq-hq-anammox-clustered.aln"
    output:
        results + "{subunit}-gtdb-refseq-hq-anammox-clustered-trimal.aln"
    conda:
        "../envs/trimal.yaml"
    shell:
        "trimal -automated1 -in {input} -out {output}"

rule fasttree:
    input:
        results + "{subunit}-gtdb-refseq-hq-anammox-clustered-trimal.aln"
    output:
        results + "{subunit}-gtdb-refseq-hq-anammox-clustered-trimal.fasttree"
    conda:
        "../envs/fasttree.yaml"
    shell:
        "fasttree -lg {input} > {output}"

rule iqtree_clustered_nxr:
    input:
        results + "{subunit}-gtdb-refseq-hq-anammox-clustered-trimal.aln"
    output:
        results + "{subunit}-gtdb-refseq-hq-anammox-clustered-trimal.treefile"
    params:
        pre = results + "{subunit}-gtdb-refseq-hq-anammox-clustered-trimal"
    threads:
        12
    conda:
        envs + "iqtree.yaml"
    shell:
        "iqtree2 -s {input} -m MFP -mset LG,WAG -B 1000 -alrt 1000 -pre {params.pre} -nt {threads} -redo"

rule add_tree_annotation_column:
    input:
        local_taxonomy = data + "local_taxonomy.tsv",
        local_seq = data + "nxr-phylogeny/{subunit}-hq-anammox.faa",
        table = results + "{subunit}-gtdb-refseq-no-anammox.tsv"
    output:
        results + "{subunit}-phylo-annotation.tsv"
    conda:
        envs + "biopython.yaml"
    shell:
        "python ../general/merge-blast-table-local-taxa.py {input.table} {input.local_taxonomy} {input.local_seq} {output}"

rule extract_nxr_a_clade_accessions:
    input:
        clusters = results + "{subunit}-gtdb-refseq-hq-anammox-clustered.faa.clstr",
        cluster_accessions = data + "nxr-phylogeny/{subunit}-clade-accessions.txt"
    output:
        results + "{subunit}-gtdb-refseq-hq-anammox-clade-accessions.txt"
    shell:
        "python extract-clustered-sequences.py {input.cluster_accessions} {input.clusters} > {output}"

rule extract_nxr_clade_outgroup_sequences:
    input:
        clade_accessions = results + "{subunit}-gtdb-refseq-hq-anammox-clade-accessions.txt",
        outgroup_accessions = data + "nxr-phylogeny/{subunit}-outgroup-accessions.txt",
        sequences = results + "{subunit}-gtdb-refseq-hq-anammox.faa"
    output:
        results + "{subunit}-gtdb-refseq-clade-outgroup.faa"
    conda:
        envs + "seqtk.yaml"
    shell:
        "cat {input.clade_accessions} {input.outgroup_accessions} | seqtk subseq {input.sequences} - > {output}"

rule align_nxr_clade_outgroup:
    input:
        results + "{subunit}-gtdb-refseq-clade-outgroup.faa"
    output:
        results + "{subunit}-gtdb-refseq-clade-outgroup.aln"
    conda:
        envs + "mafft.yaml"
    threads:
        12
    shell:
        "mafft-linsi --reorder --thread {threads} {input} > {output}"

rule trim_nxr_clade_outgroup:
    input:
        results + "{subunit}-gtdb-refseq-clade-outgroup.aln"
    output:
        results + "{subunit}-gtdb-refseq-clade-outgroup-trimal.aln"
    conda:
        envs + "trimal.yaml"
    shell:
        "trimal -automated1 -in {input} -out {output}"

rule iqtree_nxr_clade_outgroup:
    input:
        results + "{subunit}-gtdb-refseq-clade-outgroup-trimal.aln"
    output:
        results + "{subunit}-gtdb-refseq-clade-outgroup-trimal.treefile"
    params:
        pre = results + "{subunit}-gtdb-refseq-clade-outgroup-trimal"
    conda:
        envs + "iqtree.yaml"
    threads:
        12
    shell:
        "iqtree2 -s {input} -m MFP -mset LG,WAG -B 1000 -alrt 1000 -pre {params.pre} -nt {threads} -redo"

rule extract_nxr_c_gtdb:
    input:
        results + "nxr-gtdb-w-taxa.tsv.gz"
    output:
        sequences=results + "nxr_c-gtdb.faa",
        table=results + "nxr_c-gtdb.tsv"
    shell:
        """
        python ../general/extract-blast-hits.py \
            --blast {input} \
            --fasta {output.sequences} \
            --table {output.table} \
            --queries SOH03958.1 \
            --qcover 0.5 \
            --scover 0.5 \
            --anammox
        """

rule extract_nxr_d_gtdb:
    input:
        results + "nxr-gtdb-w-taxa.tsv.gz"
    output:
        sequences=results + "nxr_d-gtdb.faa",
        table=results + "nxr_d-gtdb.tsv"
    shell:
        """
        python ../general/extract-blast-hits.py \
            --blast {input} \
            --fasta {output.sequences} \
            --table {output.table} \
            --queries SOH03956.1 \
            --qcover 0.5 \
            --scover 0.5 \
            --anammox
        """

rule extract_nxr_hypo_gtdb:
    input:
        results + "nxr-gtdb-w-taxa.tsv.gz"
    output:
        sequences=results + "nxr_hypo-gtdb.faa",
        table=results + "nxr_hypo-gtdb.tsv"
    shell:
        """
        python ../general/extract-blast-hits.py \
            --blast {input} \
            --fasta {output.sequences} \
            --table {output.table} \
            --queries SOH03955.1 \
            --qcover 0.5 \
            --scover 0.5 \
            --anammox
        """

rule cluster_nxr_gtdb:
    input:
        results + "{subunit}-gtdb.faa"
    output:
        results + "{subunit}-gtdb-clustered.faa"
    conda:
        "../envs/cd-hit.yaml"
    shell:
        "cd-hit -i {input} -o {output} -c 0.9"

rule align_gtdb:
    input:
        results + "{subunit}-gtdb-clustered.faa"
    output:
        results + "{subunit}-gtdb-clustered.aln"
    conda:
        "../envs/mafft.yaml"
    threads:
        12
    shell:
        "mafft-linsi --thread {threads} --reorder {input} > {output}"

rule trim_alignment_gtdb:
    input:
        results + "{subunit}-gtdb-clustered.aln"
    output:
        results + "{subunit}-gtdb-clustered-trimal.aln"
    conda:
        "../envs/trimal.yaml"
    shell:
        "trimal -automated1 -in {input} -out {output}"

rule fasttree_gtdb:
    input:
        results + "{subunit}-gtdb-clustered-trimal.aln"
    output:
        results + "{subunit}-gtdb-clustered-trimal.fasttree"
    conda:
        "../envs/fasttree.yaml"
    shell:
        "fasttree -lg {input} > {output}"
