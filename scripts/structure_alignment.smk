"""
Workflow to perform structure alignments of the HZS proteins.
"""
results = "../processed_data/hzs-structure-alignment-no-masking"
backup_dir = "/media/argos-emiha442/emiha442/1_projects/1_221123_anammox_pathway_evo/processed_data/hzs-structure-alignment-no-masking/"

rule all:
    input:
        expand(results + "/{subunit}_alphafold_identifiers.tsv", subunit=["hzs_bc", "hzs_a"]),
        results + "/hzs_bc_structures",
        results + "/hzs_a_structures",
        results + "/hzs_bc_af_structure.pdb",
        results + "/hzs_a_af_structure.pdb",
        os.path.join(backup_dir, "alignment_summary_plot.pdf"),
        os.path.join(backup_dir, "alignment_summary.tsv")


rule find_alpha_fold_structures:
    input:
        blast_df="../processed_data/hzs-outside-anammox-no-masking/hzs.gtdb.w_taxa.tsv.gz"
    output:
        results + "/{subunit}_alphafold_identifiers.tsv"
    params:
        subunit="{subunit}"
    shell:
        "python search_structures_uniprot.py {input} {output} {params.subunit}"

rule download_hzs_a_structure:
    output:
        results + "/hzs_a_af_structure.pdb"
    shell:
        "wget -O {output} https://alphafold.ebi.ac.uk/files/AF-A0A286U438-F1-model_v4.pdb"

rule download_hzs_bc_structure:
    output:
        results + "/hzs_bc_af_structure.pdb"
    shell:
        "wget -O {output} https://alphafold.ebi.ac.uk/files/AF-A0A286U486-F1-model_v4.pdb"

checkpoint download_homolog_hzs_bc_structures:
    input:
        results + "/hzs_bc_alphafold_identifiers.tsv"
    output:
        directory(results + "/hzs_bc_structures")
    shell:
        """
        awk -F'\\t' '{{ print $2 }}' {input} | grep -v 'alphafold_url' | sed '/^$/d' | wget -i - -P {output}
        """


def get_hzs_bc_structures(wildcards):
    ck_output = checkpoints.download_homolog_hzs_bc_structures.get(**wildcards).output[0]
    homologs, = glob_wildcards(os.path.join(ck_output, "{homologs}.pdb"))
    return os.path.join(results + "/hzs_bc_structures",
                               "{homologs}.pdb")

rule structure_alignment_hzs_bc:
    input:
        ref=results + "/hzs_bc_af_structure.pdb",
        target=get_hzs_bc_structures
    output:
        results + "/hzs_bc_structure_alignment/{homologs}.struct_aln"
    conda:
        "envs/tmalign.yaml"
    shell:
        """
        TMalign {input.ref} {input.target} -a > {output}
        """

checkpoint download_homolog_hzs_a_structures:
    input:
        results + "/hzs_a_alphafold_identifiers.tsv"
    output:
        directory(results + "/hzs_a_structures")
    shell:
        """
        awk -F'\\t' '{{ print $2 }}' {input} | grep -v 'alphafold_url' | sed '/^$/d' | wget -i - -P {output}
        """

def get_hzs_a_structures(wildcards):
    ck_output = checkpoints.download_homolog_hzs_a_structures.get(**wildcards).output[0]
    homologs, = glob_wildcards(os.path.join(ck_output, "{homologs}.pdb"))
    return os.path.join(results + "/hzs_a_structures",
                               "{homologs}.pdb")

rule structure_alignment_hzs_a:
    input:
        ref=results + "/hzs_a_af_structure.pdb",
        target=get_hzs_a_structures
    output:
        results + "/hzs_a_structure_alignment/{homologs}.struct_aln"
    conda:
        "envs/tmalign.yaml"
    shell:
        """
        TMalign {input.ref} {input.target} -a > {output}
        """

def get_hzs_bc_structure_alignment(wildcards):
    ck_output = checkpoints.download_homolog_hzs_bc_structures.get(**wildcards).output[0]
    homologs, = glob_wildcards(os.path.join(ck_output, "{homologs}.pdb"))
    return expand(os.path.join(results + "/hzs_bc_structure_alignment", "{homologs}.struct_aln"), homologs=homologs)

def get_hzs_a_structure_alignment(wildcards):
    ck_output = checkpoints.download_homolog_hzs_a_structures.get(**wildcards).output[0]
    homologs, = glob_wildcards(os.path.join(ck_output, "{homologs}.pdb"))
    return expand(os.path.join(results + "/hzs_a_structure_alignment", "{homologs}.struct_aln"), homologs=homologs)

rule summarize_results:
    input:
        get_hzs_bc_structure_alignment,
        get_hzs_a_structure_alignment
    output:
        results + "/alignment_summary.tsv"
    shell:
        "python parse_structure_alignments.py {input} {output}"

rule plot_results:
    input:
        results + "/alignment_summary.tsv"
    output:
        results + "/alignment_summary_plot.pdf"
    shell:
        "python plot_structure_alignment_summary.py {input} {output}"

rule backup_results:
    input:
        results + "/alignment_summary.tsv"
    output:
        os.path.join(backup_dir, "alignment_summary.tsv")
    shell:
        "cp {input} {output}"

rule backup_figure:
    input:
        results + "/alignment_summary_plot.pdf"
    output:
        os.path.join(backup_dir, "alignment_summary_plot.pdf")
    shell:
        "cp {input} {output}"
