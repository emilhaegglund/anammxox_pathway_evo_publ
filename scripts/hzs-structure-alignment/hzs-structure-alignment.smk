"""
Workflow to perform structure alignments of the HZS proteins.
"""
results = "../../processed-data/hzs-structure-alignment/"
backup_dir = "/media/argos-emiha442/emiha442/1_projects/1_221123_anammox_pathway_evo/processed_data/hzs-structure-alignment-no-masking/"
envs = "../envs/"

rule all:
    input:
        expand(results + "{subunit}-alphafold-identifiers.tsv", subunit=["hzs_bc", "hzs_a"]),
        results + "hzs_bc-structures",
        results + "hzs_a-structures",
        results + "hzs_bc-af-structure.pdb",
        results + "hzs_a-af-structure.pdb",
        os.path.join(backup_dir, "alignment-summary.pdf"),
        os.path.join(backup_dir, "alignment-summary.tsv")


rule find_alpha_fold_structures:
    """
    Parse uniprot and collect links to AlphaFold structures for the hits
    which have structures available.
    """
    input:
        blast_df="../../processed-data/hzs-homologs/hzs.gtdb.w_taxa.tsv.gz"
    output:
        results + "/{subunit}-alphafold-identifiers.tsv"
    params:
        subunit="{subunit}"
    shell:
        "python search-structures-uniprot.py {input} {output} {params.subunit}"

rule download_hzs_a_structure:
    """
    Download Scalinduae 
    """
    output:
        results + "hzs_a-af-structure.pdb"
    shell:
        "wget -O {output} https://alphafold.ebi.ac.uk/files/AF-A0A286U438-F1-model_v4.pdb"

rule download_hzs_bc_structure:
    """
    Download Scalinduae 
    """
    output:
        results + "hzs_bc-af-structure.pdb"
    shell:
        "wget -O {output} https://alphafold.ebi.ac.uk/files/AF-A0A286U486-F1-model_v4.pdb"

checkpoint download_homolog_hzs_bc_structures:
    input:
        results + "hzs_bc-alphafold-identifiers.tsv"
    output:
        directory(results + "hzs_bc-structures")
    shell:
        """
        awk -F'\\t' '{{ print $2 }}' {input} | grep -v 'alphafold_url' | sed '/^$/d' | wget -i - -P {output}
        """

def get_hzs_bc_structures(wildcards):
    ck_output = checkpoints.download_homolog_hzs_bc_structures.get(**wildcards).output[0]
    homologs, = glob_wildcards(os.path.join(ck_output, "{homologs}.pdb"))
    return os.path.join(results + "hzs_bc-structures",
                               "{homologs}.pdb")

rule structure_alignment_hzs_bc:
    input:
        ref=results + "hzs_bc-af-structure.pdb",
        target=get_hzs_bc_structures
    output:
        results + "hzs_bc-structure-alignment/{homologs}.txt"
    conda:
        envs + "tmalign.yaml"
    shell:
        """
        TMalign {input.ref} {input.target} -a > {output}
        """

checkpoint download_homolog_hzs_a_structures:
    input:
        results + "hzs_a-alphafold-identifiers.tsv"
    output:
        directory(results + "hzs_a-structures")
    shell:
        """
        awk -F'\\t' '{{ print $2 }}' {input} | grep -v 'alphafold_url' | sed '/^$/d' | wget -i - -P {output}
        """

def get_hzs_a_structures(wildcards):
    ck_output = checkpoints.download_homolog_hzs_a_structures.get(**wildcards).output[0]
    homologs, = glob_wildcards(os.path.join(ck_output, "{homologs}.pdb"))
    return os.path.join(results + "hzs_a-structures",
                               "{homologs}.pdb")

rule structure_alignment_hzs_a:
    input:
        ref=results + "hzs_a-af-structure.pdb",
        target=get_hzs_a_structures
    output:
        results + "hzs_a-structure-alignment/{homologs}.txt"
    conda:
        envs + "tmalign.yaml"
    shell:
        """
        TMalign {input.ref} {input.target} -a > {output}
        """

def get_hzs_bc_structure_alignment(wildcards):
    ck_output = checkpoints.download_homolog_hzs_bc_structures.get(**wildcards).output[0]
    homologs, = glob_wildcards(os.path.join(ck_output, "{homologs}.pdb"))
    return expand(os.path.join(results + "hzs_bc-structure-alignment", "{homologs}.txt"), homologs=homologs)

def get_hzs_a_structure_alignment(wildcards):
    ck_output = checkpoints.download_homolog_hzs_a_structures.get(**wildcards).output[0]
    homologs, = glob_wildcards(os.path.join(ck_output, "{homologs}.pdb"))
    return expand(os.path.join(results + "hzs_a-structure-alignment", "{homologs}.txt"), homologs=homologs)

rule summarize_results:
    input:
        get_hzs_bc_structure_alignment,
        get_hzs_a_structure_alignment
    output:
        results + "alignment-summary.tsv"
    shell:
        "python parse_structure_alignments.py {input} {output}"

rule plot_results:
    input:
        results + "alignment-summary.tsv"
    output:
        results + "alignment-summary.pdf"
    shell:
        "python plot-structure-alignment-summary.py {input} {output}"

rule backup_results:
    input:
        results + "alignment_summary.tsv"
    output:
        os.path.join(backup_dir, "alignment-summary.tsv")
    shell:
        "cp {input} {output}"

rule backup_figure:
    input:
        results + "alignment-summary.pdf"
    output:
        os.path.join(backup_dir, "alignment-summary.pdf")
    shell:
        "cp {input} {output}"
