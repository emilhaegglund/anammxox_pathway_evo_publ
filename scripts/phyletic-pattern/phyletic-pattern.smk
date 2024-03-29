results = "../../processed-data/phyletic-pattern/"
backup_dir = "/media/argos-emiha442/emiha442/1_projects/1_221123_anammox_pathway_evo/processed-data/phyletic-pattern/"
envs = "../envs/"
data = "../../data/"
# Get proteomes to use
PVC_ACCESSIONS = []
with open(data + "rpoB-phylogeny/pvc-accessions.txt", "r") as f:
    for line in f:
        PVC_ACCESSIONS.append(line.strip())

ANAMMOX_ACCESSIONS, = glob_wildcards(data + "anammox-proteomes/{accession}_protein.faa")

rule all:
    input:
        results + "anammox-pathway-kust-blastp-pvc.tsv",
        results + "pvc-anammox.dmnd",
        backup_dir + "phyletic-pattern.pdf",
        backup_dir + "phyletic-pattern-supplementary.pdf",
        results + "anammox-pathway-kust.faa"

rule makedb:
    input:
        pvc=expand("../../processed-data/rpoB-phylogeny/assembly-data/ncbi_dataset/data/{accession}/protein.faa", accession=PVC_ACCESSIONS),
        anammox=expand(data + "anammox-proteomes/{accession}_protein.faa", accession=ANAMMOX_ACCESSIONS)
    output:
        results + "pvc-anammox.dmnd"
    params:
        prefix = results + "pvc-anammox"
    conda:
        envs + "diamond.yaml"
    shell:
        "cat {input.pvc} {input.anammox} | diamond makedb --db {params.prefix}"

rule download_queries:
    input:
        data + "phyletic-pattern/anammox-pathway-proteins.txt"
    output:
        results + "anammox-pathway-kust.faa"
    conda:
        envs + "entrez.yaml"
    shell:
       "grep -v '#' {input} | awk -F'\\t' '{{ print $3 }}' | sed 1d | efetch -db protein -format fasta > {output}"

rule diamond:
    input:
        query = results + "anammox-pathway-kust.faa",
        db = results + "pvc-anammox.dmnd"
    output:
        results + "anammox-pathway-kust-blastp-pvc.tsv"
    conda:
        envs + "diamond.yaml"
    threads:
        24
    shell:
        """
        diamond blastp --db {input.db} \
            --query {input.query} \
            --out {output} \
            --outfmt 6 \
            --ultra-sensitive \
            --evalue 1e-6 \
            --max-target-seqs 500 \
            --threads {threads}
        """

rule map_protein_to_accession:
    input:
        pvc_proteomes = "../../processed-data/rpoB-phylogeny/protein-mappings.tsv",
        anammox_proteomes = data + "anammox-proteomes"
    output:
        results + "proteins-genome-accessions.tsv"
    conda:
        envs + "biopython.yaml"
    shell:
        "python map-protein-to-genome.py {input.pvc_proteomes} {input.anammox_proteomes} {output}"

rule plot_supplementary_figure:
    input:
        sort_order = data + "phyletic-pattern/anammox-extended-rpoB-order.txt",
        queries = results + "anammox-pathway-kust.faa",
        protein_accession = results + "proteins-genome-accessions.tsv",
        blast = results + "anammox-pathway-kust-blastp-pvc.tsv",
        pathway_proteins = data + "phyletic-pattern/anammox-pathway-proteins.txt",
        genome_metadata = "../../processed-data/rpoB-phylogeny/rpoB-tree-annotation.tsv",
    output:
        backup_dir + "phyletic-pattern-supplementary.pdf"
    params:
        width=8,
        height=24
    conda:
        envs + "biopython.yaml"
    shell:
        """
        python plot-presence-absence.py \
            {input.sort_order} \
            {input.queries} \
            {input.protein_accession} \
            {input.blast} \
            {input.pathway_proteins} \
            {input.genome_metadata} \
            {output} \
            {params.width} \
            {params.height}
        """

rule plot_figure:
    input:
        sort_order = data + "phyletic-pattern/anammox-hq-rpoB-order.txt",
        queries = results + "anammox-pathway-kust.faa",
        protein_accession = results + "proteins-genome-accessions.tsv",
        blast = results + "anammox-pathway-kust-blastp-pvc.tsv",
        pathway_proteins=data + "phyletic-pattern/anammox-pathway-proteins.txt",
        genome_metadata="../../processed-data/rpoB-phylogeny/rpoB-tree-annotation.tsv"
    output:
        backup_dir + "phyletic-pattern.pdf"
    params:
        width=8,
        height=8
    conda:
        envs + "biopython.yaml"
    shell:
        """
        python plot-presence-absence.py \
            {input.sort_order} \
            {input.queries} \
            {input.protein_accession} \
            {input.blast} \
            {input.pathway_proteins} \
            {input.genome_metadata} \
            {output} \
            {params.width} \
            {params.height}
        """

#rule backup_figure:
#    input:
#        results + "phyletic-pattern.pdf"
#    output:
#        backup_dir + "phyletic-pattern.pdf"
#    shell:
#        "cp {input} {output}"

#rule backup_supplementary_figure:
#    input:
#        results + "phyletic-pattern-supplementary.pdf"
#    output:
#        backup_dir + "phyletic-pattern-supplementary.pdf"
#    shell:
#        "cp {input} {output}"
