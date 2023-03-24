
rule all:
    input:
        "../processed_data/hzs-recomb-dnds/hzs-a-codeml/hzs-a.codeml.filtered.txt"

rule codeml:
    input:
        "../data/hzs-a.codeml.ctl"
    output:
        "../processed_data/hzs-recomb-dnds/hzs-a-codeml/hzs-a.codeml.txt"
    params:
        outdir="../processed_data/hzs-recomb-dnds/hzs-a-codeml/"
    conda:
        "envs/paml.yaml"
    shell:
        """
        cp {input} {params.outdir}/codeml.ctl;
        cd {params.outdir};
        codeml;
        """
rule parse_codeml:
    input:
        "../processed_data/hzs-recomb-dnds/hzs-a-codeml/hzs-a.codeml.txt"
    output:
        "../processed_data/hzs-recomb-dnds/hzs-a-codeml/hzs-a.codeml.filtered.txt"
    shell:
        "python parse_codeml_output.py {input} > {output}"
