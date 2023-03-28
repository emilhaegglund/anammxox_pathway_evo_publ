backup_dir = "/media/argos-emiha442/emiha442/1_projects/1_221123_anammox_pathway_evo/processed_data/nxr-phylogeny/"
results = "../../processed_data/nxr-phylogeny/"
envs = "../envs/"
rule all:
    input:
        expand(results + "{subunit}.gtdb.refseq.hits.no_anammox.faa", subunit=["nxr-a", "nxr-b"])

rule download_queries:
    output:
        results + "nxr.query.faa"
    conda:
        envs + "entrez.yaml"
    shell:
        "efetch -id SOH03954.1,SOH03957.1,SOH03958.1 -db protein -format fasta > {output}"

rule blast_gtdb:
    input:
        results + "nxr.query.faa"
    output:
        results + "nxr.gtdb.tsv.gz"
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
        blast=results + "nxr.gtdb.tsv.gz",
        prot="../../data/gtdb_representatives_database/prot2accession.tsv.gz"
    output:
        temp(results + "nxr.gtdb.w_accessions.tsv.gz")
    shell:
        "python ../hzs-outside-anammox-merge-blast-accessions.py {input.blast} {input.prot} {output}"

rule add_taxa:
    input:
        blast=results + "nxr.gtdb.w_accessions.tsv.gz",
        gtdb="../../data/gtdb_representatives_database/gtdb_representatives.sample_gtdb.metadata.tsv",
    output:
        results + "nxr.gtdb.w_taxa.tsv.gz"
    shell:
        "python ../hzs-outside-anammox-merge-blast-taxa.py {input.blast} {input.gtdb} {output}"

rule extract_nxr_a_gtdb_refseq:
    input:
        results + "nxr.gtdb.w_taxa.tsv.gz"
    output:
        results + "nxr-a.gtdb.refseq.hits.no_anammox.faa"
    shell:
        """python ../extract-blast-hits.py --blast {input} \
            --fasta {output} \
            --queries SOH03954.1 \
            --qcover 0.5 \
            --scover 0.5 \
            --source refseq
        """

rule extract_nxr_b_gtdb_refseq:
    input:
        results + "nxr.gtdb.w_taxa.tsv.gz"
    output:
        results + "nxr-b.gtdb.refseq.hits.no_anammox.faa"
    shell:
        """python ../extract-blast-hits.py --blast {input} \
            --fasta {output} \
            --queries SOH03957.1 \
            --qcover 0.5 \
            --scover 0.5 \
            --source refseq
        """