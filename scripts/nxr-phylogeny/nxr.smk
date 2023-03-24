backup_dir = "/media/argos-emiha442/emiha442/1_projects/1_221123_anammox_pathway_evo/processed_data/nxr-phylogeny/"
results = "../../processed_data/nxr-phylogeny/"
envs = "../envs/"
rule all:
    input:
        results + "nxr.gtdb.tsv.gz"

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