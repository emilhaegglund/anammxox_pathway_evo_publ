"""
Workflow to map the proteins of Ca. Kuenenia stuttgartiensis MBR1 to other
isolates of Ca. Kuenenia stuttgartiensis.
"""
rule all:
    input:
        expand("../processed_data/kuenenia_mappings/MBR1_{proteome}.blastp.tsv",
        proteome=["kust", "CSTR1_RefSeq"]),
        "../processed_data/kuenenia_mappings/MBR1_kust.locus_tags.blastp.tsv"

rule download_CSTR1_RefSeq:
    output:
        "../processed_data/kuenenia_mappings/proteomes/CSTR1_RefSeq.faa"
    conda:
        "envs/ncbi-datasets.yaml"
    shell:
        """
        datasets download genome accession GCF_011066545.1 \
            --include protein \
            --filename ../processed_data/kuenenia_mappings/proteomes/GCF_011066545.1.zip;
        unzip -p ../processed_data/kuenenia_mappings/proteomes/GCF_011066545.1.zip \
            ncbi_dataset/data/GCF_011066545.1/protein.faa > {output};
        """


rule download_kust_proteome:
    """
    Download the proteins from the kust assembly.
    """
    output:
        "../processed_data/kuenenia_mappings/proteomes/kust.faa",
    conda:
       "envs/entrez.yaml"
    threads:
        1
    shell:
        """
        esearch -db nuccore -query CT030148.2 | elink -target protein |\
            efetch -format fasta >> {output};
        esearch -db nuccore -query CT573074.1 | elink -target protein |\
            efetch -format fasta >> {output};
        esearch -db nuccore -query CT573073.1 | elink -target protein |\
            efetch -format fasta >> {output};
        esearch -db nuccore -query CT573072.1 | elink -target protein |\
            efetch -format fasta >> {output};
        esearch -db nuccore -query CT573071.1 | elink -target protein |\
            efetch -format fasta >> {output};
        """

rule create_diamond_db:
    input:
        "../processed_data/kuenenia_mappings/proteomes/{proteome}.faa"
    output:
        "../processed_data/kuenenia_mappings/{proteome}.dmnd"
    shell:
        "cat {input} | diamond makedb --db {output}"

rule diamond:
    input:
        query="../data/proteomes/GCA_900232105.1_Kuenenia_stuttgartiensis_MBR1_protein.faa",
        db="../processed_data/kuenenia_mappings/{proteome}.dmnd"
    output:
        "../processed_data/kuenenia_mappings/MBR1_{proteome}.blastp.tsv"
    conda:
        "envs/diamond.yaml"
    shell:
        """
        diamond blastp --db {input.db} \
        --query {input.query} \
        --out {output} \
        --outfmt 6 \
        --query-cover 70 \
        --subject-cover 70 \
        --max-target-seqs 1 \
        --sensitive \
        --threads {threads}
        """

rule download_kust_genbank:
    """
    Download the genbank to map protein accessions to locus tags for kust.
    """
    output:
        "../processed_data/kuenenia_mappings/genbanks/kust.gbk",
    conda:
       "envs/entrez.yaml"
    threads:
        1
    shell:
        """
        esearch -db nuccore -query CT030148.2 | \
            efetch -format genbank >> {output};
        esearch -db nuccore -query CT573074.1 | \
            efetch -format genbank >> {output};
        esearch -db nuccore -query CT573073.1 | \
            efetch -format genbank >> {output};
        esearch -db nuccore -query CT573072.1 | \
            efetch -format genbank >> {output};
        esearch -db nuccore -query CT573071.1 | \
            efetch -format genbank >> {output};
        """

rule add_kust_locus_tags:
    input:
        blast="../processed_data/kuenenia_mappings/MBR1_kust.blastp.tsv",
        genbank="../processed_data/kuenenia_mappings/genbanks/kust.gbk"
    output:
        "../processed_data/kuenenia_mappings/MBR1_kust.locus_tags.blastp.tsv"
    conda:
        "envs/biopython.yaml"
    shell:
        "python map-kust-protein-locus-tags.py {input.blast} {input.genbank} {output}"
