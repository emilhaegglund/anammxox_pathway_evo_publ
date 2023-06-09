"""
Workflow to produce the data needed for the HZS-operon overview figure.
"""

backup_dir = "/media/argos-emiha442/emiha442/1_projects/1_221123_anammox_pathway_evo/processed-data/hzs-operon-overview/"
results = "../../processed-data/hzs-operon-overview/"
processed_data = "../../processed-data/"
data = "../../data/"
envs = "../envs/"


rule all:
    input:
        results + "hzs-orthogroups.tsv",
        results + "hzs-operon.tsv",
        results + "hzs-operon-proteins.blastp.tsv",
        backup_dir + "hzs-operon.tsv",
        backup_dir + "hzs-orthogroups.tsv",
        backup_dir + "hzs-operon-blast.tsv"

rule extract_data:
    input:
        orthogroups = processed_data + "orthofinder/Results_Jun08/Orthogroups/Orthogroups.txt",
        genbank_dir = data + "anammox-genbanks"
    output:
        orthogroups = results + "hzs-orthogroups.tsv",
        gene_positions = results + "hzs-operon.tsv"
    conda:
        envs + "biopython.yaml"
    shell:
        "python hzs-operon-overview.py {input.orthogroups} {input.genbank_dir} {output.orthogroups} {output.gene_positions} "

rule concat_anammox_proteins:
    input:
        data + "anammox-proteomes/"
    output:
        temp(results + "all_anammox.faa")
    shell:
        "cat {input}/*.faa > {output}"

rule extract_proteins:
    input:
        operon_data=results + "hzs-operon.tsv",
        all_proteins=results + "all_anammox.faa"
    output:
        results + "hzs-operon-proteins.faa"
    conda:
        envs + "seqtk.yaml"
    shell:
        "awk -F'\\t' '{{ print $2 }}' {input.operon_data} | seqtk subseq {input.all_proteins} - > {output}"

rule diamond_database:
    input:
        results + "hzs-operon-proteins.faa"
    output:
        results + "hzs-operon-proteins.dmnd"
    conda:
        envs + "diamond.yaml"
    shell:
        "diamond makedb --db {output} --in {input}"

rule blast_all_against_all:
    input:
        db = results + "hzs-operon-proteins.dmnd",
        queries = results + "hzs-operon-proteins.faa"
    output:
        results + "hzs-operon-proteins.blastp.tsv"
    conda:
        envs + "diamond.yaml"
    threads:
        24
    shell:
        """
        diamond blastp --db {input.db} \
                        --query {input.queries} \
                        --out {output} \
                        --outfmt 6 \
                        --max-target-seqs 50 \
                        --threads {threads} \
                        --ultra-sensitive
        """

rule prepare_blast_results:
    input:
        blast = "../prochzs-operon-proteins.blastp.tsv",
        operon_data="../processed_data/hzs-operon-overview/hzs-operon.tsv"
    output:
        "../processed_data/hzs-operon-overview/hzs-operon-blast.tsv"
    shell:
        "python hzs-operon-overview-prepare-blast.py {input.blast} {input.operon_data} {output}"

rule backup_orthogroup:
    input:
        results + "hzs-orthogroups.tsv",
    output:
        backup_dir + "hzs-orthogroups.tsv"
    shell:
        "cp {input} {output}"

rule backup_gene_postions:
    input:
        results + "hzs-operon.tsv"
    output:
        backup_dir + "/hzs-operon.tsv"
    shell:
        "cp {input} {output}"

rule backup_blast_data:
    input:
        results + "hzs-operon-blast.tsv"
    output:
        backup_dir + "hzs-operon-blast.tsv"
    shell:
        "cp {input} {output}"
