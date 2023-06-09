backup_dir = "/media/argos-emiha442/emiha442/1_projects/1_221123_anammox_pathway_evo/processed-data/hzs-operon-phylogenies/"
results = "../../processed-data/hzs-operon-phylogenies/"
processed_data = "../../processed-data/"
envs = "../envs/"

orthogroups=["OG0000182",
             "OG0000124",
             "OG0000271",
             ]

rule all:
    input:
        expand(results + "phylogenies/{orthogroup}.operon.treefile", orthogroup=orthogroups),
        expand(backup_dir + "phylogenies/{orthogroup}.operon.treefile", orthogroup=orthogroups),
        expand(backup_dir + "alignments/{orthogroup}.aln", orthogroup=orthogroups)

rule concat_anammox_proteins:
    input:
        "../../data/anammox-proteomes/"
    output:
        temp(results + "all_anammox.faa")
    shell:
        "cat {input}/*.faa > {output}"

rule select_sequences:
    input:
        operon_orthogroups = processed_data + "hzs-operon-overview/hzs-orthogroups.tsv",
        all_anammox = results + "all_anammox.faa"
    output:
        results + "sequences/{orthogroup}.fa"
    params:
        orthogroup="{orthogroup}"
    conda:
        envs + "biopython.yaml"
    shell:
        "python hzs-operon-phylogenies-sequences.py {params.orthogroup} {input.operon_orthogroups} {input.all_anammox} {output}"

rule align:
    input:
        results + "sequences/{orthogroup}.fa"
    output:
        results + "alignments/{orthogroup}.aln"
    conda:
        envs + "mafft.yaml"
    threads:
        12
    shell:
        "mafft-linsi --thread {threads} {input} > {output}"

rule trim_alignment:
    input:
        results + "alignments/{orthogroup}.aln"
    output:
        results + "trimmed_alignments/{orthogroup}.trimal.aln"
    conda:
        envs + "trimal.yaml"
    shell:
        "trimal -automated1 -in {input} -out {output}"

rule iqtree:
    input:
        results + "trimmed_alignments/{orthogroup}.trimal.aln"
    output:
        results + "phylogenies/{orthogroup}.operon.treefile"
    params:
        prefix = results + "phylogenies/{orthogroup}.operon"
    conda:
        envs + "iqtree.yaml"
    threads:
        12
    shell:
        "iqtree2 -s {input} -m MFP -mset LG,WAG -B 1000 -pre {params.prefix} -T {threads} -redo"

rule backup_phylogenies:
    input:
        results + "phylogenies/{orthogroup}.operon.treefile"
    output:
        backup_dir + "phylogenies/{orthogroup}.operon.treefile"
    shell:
        "cp {input} {output}"

rule backup_alignments:
    input:
        results + "alignments/{orthogroup}.aln"
    output:
        backup_dir + "alignments/{orthogroup}.aln"
    shell:
        "cp {input} {output}"
