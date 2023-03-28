"""
Workflow to perform Recombination and dN/dS analysis of the HZS
"""
results = "../../processed-data/hzs-dnds/"
processed_data = "../../processed-data/"
data = "../../data/"
backup_dir = "/media/argos-emiha442/emiha442/1_projects/1_221123_anammox_pathway_evo/processed_data/hzs-dnds"

HZS_SUBUNITS = ["hzs-a", "hzs-b", "hzs-c"]

rule all:
    input:
        expand(backup_dir + "/hzs-recomb-dnds/{subunit}.nucleotide-identity.tsv", subunit=HZS_SUBUNITS),
        expand(backup_dir + "/hzs-recomb-dnds/{subunit}-nucleotide-identity-heatmap.pdf", subunit=HZS_SUBUNITS),
        expand(backup_dir + "/hzs-recomb-dnds/{subunit}-codeml/{subunit}.codeml.txt", subunit=HZS_SUBUNITS),
        expand(backup_dir + "/hzs-recomb-dnds/{subunit}-codeml/{subunit}.codeml.filtered.txt", subunit=HZS_SUBUNITS),
        expand(backup_dir + "/hzs-recomb-dnds/{subunit}-dnds-heatmap.pdf", subunit=HZS_SUBUNITS),
        backup_dir + "/hzs-recomb-dnds/hzs-dnds-violin.svg",

# Analysis of HZS Alpha subunit
rule extract_nucleotide_sequences_hzs_a:
    input:
        sequences = processed_data + "hzs-operon-phylogenies/sequences/OG0000255.fa",
        genes = data + "anammox-genes/"
    output:
        results + "hzs-a.fna"
    conda:
        envs + "biopython.yaml"
    shell:
        "python get-nucleotides.py {input.sequences} {input.genes} {output}"

rule copy_hzs_a_alignment:
    input:
        processed_data + "hzs-operon-phylogenies/alignments/OG0000255.aln"
    output:
        results + "hzs-a.aln"
    shell:
        "cp {input} {output}"

# For the analysis of the B and C subunits it's a bit more complicated since Ca.
# Scalindua japonica has this subunit fused. In the phylogenies for these subunits
# this sequence is only present in the orthogroup representing the C subunit so I manually
# added it to the orthogroup representing the B subunit and then I let the alignment trimming
# take care of the extension. For the CodeML-analysis I don't want to perform trimming since
# the pal2nal script can not handle this.
# Instead I cut the protein from Scalindua japonica so that it fits with each subunit, I also
# have to do this for the nucleotide sequence. Then add the split version of the protein to
# orthogroups.
rule cut_hzs_bc_protein:
    input:
       data + "anammox-proteomes/GCA_002443295.1_ASM244329v1_protein.faa"
    output:
        hzs_b = results + "scalindua-hzs_b.faa",
        hzs_c = results + "scalindua-hzs_c.faa"
    conda:
        envs + "biopython.yaml"
    shell:
        "python split-hzs-bc-protein.py {input} {output.hzs_b} {output.hzs_c}"

rule cut_hzs_bc_nucleotide:
    input:
       "../data/anammox-genes/GCA_002443295.1_ASM244329v1_cds_from_genomic.fna"
    output:
        hzs_b = "../processed_data/hzs-recomb-dnds/scalindua-hzs-b.fna",
        hzs_c = "../processed_data/hzs-recomb-dnds/scalindua-hzs-c.fna"
    conda:
        envs + "biopython.yaml"
    shell:
        "python split-hzs-bc-nucleotide.py {input} {output.hzs_b} {output.hzs_c}"

# Prepare HZS-B sequences for CodeML
rule remove_scalindua_hzs_b_orthogroup:
    input:
        processed_data + "hzs-operon-phylogenies/sequences/OG0000174.fa"
    output:
        results + "hzs_b-wo-scalindua.faa"
    conda:
        envs + "biopython.yaml"
    shell:
        "python remove-scalindua.py {input} {output}"


rule merge_hzs_b_proteins:
    input:
        orthogroup = results + "hzs_b-wo-scalindua.faa",
        scalindua_hzs_b = results + "scalindua-hzs_b.faa"
    output:
        results + "hzs_b.faa"
    shell:
        "cat {input.orthogroup} {input.scalindua_hzs_b} > {output}"

rule align_hzs_b:
    input:
        results + "hzs_b.faa"
    output:
        results + "hzs_b.aln"
    conda:
        envs + "mafft.yaml"
    shell:
        "mafft-linsi {input} > {output}"

rule extract_nucleotide_sequences_hzs_b:
    input:
        sequences = "../processed_data/hzs-recomb-dnds/hzs-b.wo_scalindua.faa",
        genes = data + "anammox-genes"
    output:
        results + "hzs_b-wo-scalindua.fna"
    conda:
        envs + "biopython.yaml"
    shell:
        "python get-nucleotides.py {input.sequences} {input.genes} {output}"

rule merge_hzs_b_nucleotides:
    input:
        orthogroup = results + "hzs_b-wo-scalindua.fna",
        scalindua_hzs_b = results + "scalindua-hzs_b.fna"
    output:
        results + "hzs_b.fna"
    shell:
        "cat {input.orthogroup} {input.scalindua_hzs_b} > {output}"

# Prepare HZS-C sequences for CodeML
rule remove_scalindua_hzs_c_orthogroup:
    input:
        processed_data = "hzs-operon-phylogenies/sequences/OG0000132.fa"
    output:
        results + "hzs_c-wo-scalindua.faa"
    conda:
        envs + "biopython.yaml"
    shell:
        "python remove-scalindua.py {input} {output}"

rule merge_hzs_c_proteins:
    input:
        orthogroup = results + "hzs_c-wo-scalindua.faa",
        scalindua_hzs_c = results + "scalindua-hzs_c.faa"
    output:
        results + "hzs_c.faa"
    shell:
        "cat {input.orthogroup} {input.scalindua_hzs_c} > {output}"

rule align_hzs_c:
    input:
        results + "hzs_c.faa"
    output:
        results + "hzs_c.aln"
    conda:
        envs + "mafft.yaml"
    shell:
        "mafft-linsi {input} > {output}"

rule extract_nucleotide_sequences_hzs_c:
    input:
        sequences = results + "hzs_c-wo-scalindua.faa",
        genes = data + "anammox-genes"
    output:
        results + "hzs_c-wo-scalindua.fna"
    conda:
        envs + "biopython.yaml"
    shell:
        "python get-nucleotides.py {input.sequences} {input.genes} {output}"

rule merge_hzs_c_nucleotides:
    input:
        orthogroup = results + "hzs_c-wo-scalindua.fna",
        scalindua_hzs_c = results + "scalindua-hzs_c.fna"
    output:
        results + "hzs_c.fna"
    shell:
        "cat {input.orthogroup} {input.scalindua_hzs_c} > {output}"

# Common to all subunits
rule backtranslate_aa_alignment:
    """
    Convert the AA-alignment to codon-alignment
    """
    input:
        nucleotides = results + "{subunit}.fna",
        protein_alignment = results + "{subunit}.aln"
    output:
        results + "{subunit}-codon.aln"
    conda:
        envs + "pal2nal.yaml"
    shell:
        "pal2nal.pl {input.protein_alignment} {input.nucleotides} -output fasta > {output}"

rule codeml:
    input:
        control_file = data + "hzs-recomb-dnds/{subunit}.codeml.ctl",
        codon_alignment = results + "{subunit}-codon.aln"
    output:
        results + "{subunit}-codeml/{subunit}.codeml.txt"
    params:
        outdir = results + "{subunit}-codeml"
    conda:
        envs + "paml.yaml"
    shell:
        """
        mkdir -p {params.outdir};
        cp {input.control_file} {params.outdir}/codeml.ctl;
        cd {params.outdir};
        codeml;
        """

rule parse_codeml:
    input:
        results + "{subunit}-codeml/{subunit}.codeml.txt"
    output:
        results + "{subunit}-codeml/{subunit}.codeml.filtered.txt"
    shell:
        "python parse-codeml-output.py {input} > {output}"

rule codeml_violin_plot:
    input:
        hzs_b = results + "hzs-b-codeml/hzs-b.codeml.filtered.txt",
        hzs_c = results + "hzs-c-codeml/hzs-c.codeml.filtered.txt",
        hzs_a = results + "hzs-a-codeml/hzs-a.codeml.filtered.txt",
        #hdh="../processed_data/hdh-dnds/hdh-codeml/hdh.codeml.filtered.txt",
    output:
        rersults + "hzs-dnds-violin.svg"
    shell:
        "python codeml-violin-plot.py {input.hzs_b} {input.hzs_c} {input.hzs_a} {output}"

rule codeml_heatmap:
    input:
        dnds_data = results + "{subunit}-codeml/{subunit}.codeml.filtered.txt",
        proteomes = data + "anammox-proteomes/"
        genome_metadata = data + "anammox_rpoB_subset_names.221206.tsv"
    output:
        results + "{subunit}-dnds-heatmap.pdf"
    conda:
        envs + "biopython.yaml"
    shell:
        "python codeml-heatmap.py {input.dnds_data} {input.proteomes} {input.genome_metadata} {output}"

rule calculate_nucleotide_identity:
    """
    Calculate the nucleotide identity for the pairs.
    """
    input:
        results + "hzs-recomb-dnds/{subunit}-codon.aln"
    output:
        results + "{subunit}-nucleotide-identity.tsv"
    conda:
        envs + "biopython.yaml"
    shell:
        "python nucleotide-identity.py {input} {output}"

rule nucleotide_identiy_heatmap:
    input:
        ni_data = results + "{subunit}-nucleotide-identity.tsv",
        proteomes = data + "anammox-proteomes/",
        genome_metadata = data + "anammox_rpoB_subset_names.221206.tsv"
    output:
        results + "{subunit}-nucleotide-identity-heatmap.pdf"
    conda:
        envs + "biopython.yaml"
    shell:
        "python nucleotide-identity-heatmap.py {input.ni_data} {input.proteomes} {input.genome_metadata} {output}"

# Backup results
rule backup_dnds_results:
    input:
        results + "{subunit}-codeml/{subunit}.codeml.txt"
    output:
        backup_dir + "{subunit}-codeml/{subunit}.codeml.txt"
    shell:
        "cp {input} {output}"

rule backup_dnds_filtered_values:
    input:
        results + "{subunit}-codeml/{subunit}.codeml.filtered.txt"
    output:
        backup_dir + "{subunit}-codeml/{subunit}.codeml.filtered.txt"
    shell:
        "cp {input} {output}"

rule backup_dnds_violin_plot:
    input:
        results + "hzs-dnds-violin.svg"
    output:
        backup_dir + "hzs-dnds-violin.svg"
    shell:
        "cp {input} {output}"

rule backup_dnds_heatmaps:
    input:
        results + "{subunit}-dnds-heatmap.pdf"
    output:
        backup_dir + "{subunit}-dnds-heatmap.pdf"
    shell:
        "cp {input} {output}"

rule backup_nucleotide_identity:
    input:
        results + "{subunit}-nucleotide-identity.tsv"
    output:
        backup_dir + "{subunit}-nucleotide-identity.tsv"
    shell:
        "cp {input} {output}"

rule backup_nucleotide_identity_heatmaps:
    input:
        results + "{subunit}-nucleotide-identity-heatmap.pdf"
    output:
        backup_dir + "{subunit}-nucleotide-identity-heatmap.pdf"
    shell:
        "cp {input} {output}"
