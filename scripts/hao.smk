backup_dir = "/media/argos-emiha442/emiha442/1_projects/1_221123_anammox_pathway_evo/processed_data/hao-phylogeny/"
results = "../processed_data/hao-phylogeny/"
HQ_GENOMES = ["GCA_000949635.1_ASM94963v1",
    "GCA_002443295.1_ASM244329v1",
    "GCA_007860005.1_ASM786000v1",
    "GCA_013112645.1_ASM1311264v1",
    "GCA_017347445.1_ASM1734744v1",
    "GCA_021650895.1_ASM2165089v1",
    "GCA_021650915.1_ASM2165091v1",
    "GCA_900232105.1_Kuenenia_stuttgartiensis_MBR1"]

rule all:
    input:
        results + "hao.gtdb.w_taxa.tsv.gz",
        results + "hao.hq.faa",
        results + "hao.gtdb.refseq.hq_anammox.faa",
        results + "hao.blast.annotation.tsv",
        results + "hao.gtdb.refseq.hq_anammox.trimal.fasttree",
        results + "hao.gtdb.refseq.hq_anammox.remove_long_branch_1.trimal.fasttree",
        results + "hao.gtdb.refseq.hq_anammox.remove_long_branch_1.trimal.treefile",
        #results + "hao.gtdb.refseq.hq_anammox.ougtroup.trimal.fasttree",
        #results + "hao.gtdb.refseq.hq_anammox.ougtroup.trimal.treefile",
        #results + "hao-outgroup.aln",
        #results + "hao-outgroup.anchors.faa",
        #results + "hao.gtdb.refseq.hq_anammox.remove_long_branch_1.anchors.faa",
        #results + "hao-anchors.wo_x.fasttree",
        #results + "hao-anchors.aln"
        results + "onr-outgroup-query.gtdb.blastp.tsv.gz",
        results + "onr-outgroup-query.refseq.faa",
        results + "hao.onr.trimal.treefile"
#        "../data/hq.dmnd"

rule download_queries:
    output:
        results + "query.faa"
    conda:
        "envs/entrez.yaml"
    shell:
        "efetch -id SOH04660.1,SOH06265.1,SOH05916.1,SOH05518.1,SOH05157.1,SOH04857.1,SOH03722.1,SOH04263.1,SOH05871.1,SOH05896.1 -db protein -format fasta > {output}"

rule blast_gtdb:
    input:
        results + "query.faa"
    output:
        results + "hao.gtdb.tsv.gz"
    params:
        db="../data/gtdb_representatives_database/gtdb_representatives.dmnd"
    conda:
        "envs/diamond.yaml"
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
        blast=results + "hao.gtdb.tsv.gz",
        prot="../data/gtdb_representatives_database/prot2accession.tsv.gz"
    output:
        results + "hao.gtdb.w_accessions.tsv.gz"
    shell:
        "python hzs-outside-anammox-merge-blast-accessions.py {input.blast} {input.prot} {output}"

rule add_taxa:
    input:
        blast=results + "hao.gtdb.w_accessions.tsv.gz",
        gtdb="../data/gtdb_representatives_database/gtdb_representatives.sample_gtdb.metadata.tsv",
    output:
        results + "hao.gtdb.w_taxa.tsv.gz"
    shell:
        "python hzs-outside-anammox-merge-blast-taxa.py {input.blast} {input.gtdb} {output}"

rule hq_database:
    input:
        expand("../data/proteomes/{hq}_protein.faa", hq=HQ_GENOMES)
    output:
        "../data/hq.dmnd"
    conda:
        "./envs/diamond.yaml"
    shell:
        "cat {input} | diamond makedb --db {output} --in -"

rule blast_hq_dataset:
    input:
        queries=results + "query.faa",
        db="../data/hq.dmnd"
    output:
        results + "hao.hq.blastp.tsv.gz"
    threads:
        12
    conda:
        "./envs/diamond.yaml"
    shell:
        """
        diamond blastp \
            --db {input.db} \
            --query {input.queries} \
            --subject-cover 50 \
            --query-cover 50 \
            --max-target-seqs 2000 \
            --outfmt 6 qseqid sseqid pident positive mismatch gapopen qlen slen length qstart qend sstart send evalue bitscore stitle full_sseq \
            --very-sensitive \
            --evalue 1e-6 \
            --min-score 80 \
            --out {output} \
            --compress 1 \
            --threads {threads} \
        """

rule extract_hq_hao_hits:
    input:
        results + "hao.hq.blastp.tsv.gz"
    output:
        results + "hao.hq.faa"
    shell:
        "python hao_extract_hq_hao.py {input} > {output}"

rule extract_hao_gtdb_refseq:
    input:
        results + "hao.gtdb.w_taxa.tsv.gz"
    output:
        results + "hao.gtdb.refseq.hits.no_anammox.faa"
    shell:
        """python extract-blast-hits.py --blast {input} \
            --fasta {output} \
            --qcover 0.5 \
            --scover 0.5 \
            --source refseq
        """

rule merge_anammox_hq_gtdb_refseq:
    input:
        anammox_hq=results + "hao.hq.faa",
        no_anammox=results + "hao.gtdb.refseq.hits.no_anammox.faa"
    output:
        results + "hao.gtdb.refseq.hq_anammox.faa"
    shell:
        "cat {input.anammox_hq} {input.no_anammox} > {output}"

rule align:
    input:
        results + "hao.gtdb.refseq.hq_anammox.faa"
    output:
        results + "hao.gtdb.refseq.hq_anammox.aln"
    conda:
        "./envs/mafft.yaml"
    threads:
        12
    shell:
        "mafft-linsi --thread {threads} {input} > {output}"

rule trim_alignment:
    input:
        results + "hao.gtdb.refseq.hq_anammox.aln"
    output:
        results + "hao.gtdb.refseq.hq_anammox.trimal.aln"
    conda:
        "./envs/trimal.yaml"
    shell:
        "trimal -automated1 -in {input} -out {output}"

rule fasttree:
    input:
        results + "hao.gtdb.refseq.hq_anammox.trimal.aln"
    output:
        results + "hao.gtdb.refseq.hq_anammox.trimal.fasttree"
    conda:
        "./envs/fasttree.yaml"
    shell:
        "fasttree -lg {input} > {output}"

rule remove_long_branch_1:
    input:
        sequences=results + "hao.gtdb.refseq.hq_anammox.faa",
        long_branches="../data/hao_phylogeny/remove_long_branches.txt"
    output:
        results + "hao.gtdb.refseq.hq_anammox.remove_long_branch_1.faa",
    conda:
        "./envs/biopython.yaml"
    shell:
        "python rpoB-phylogeny-remove-taxa.py {input.sequences} {input.long_branches} {output}"

rule align_after_long_branch_remove_1:
    input:
        results + "hao.gtdb.refseq.hq_anammox.remove_long_branch_1.faa",
    output:
        results + "hao.gtdb.refseq.hq_anammox.remove_long_branch_1.aln",
    conda:
        "./envs/mafft.yaml"
    threads:
        12
    shell:
        "mafft-linsi --thread {threads} {input} > {output}"

rule trim_after_long_branch_remove_1:
    input:
        results + "hao.gtdb.refseq.hq_anammox.remove_long_branch_1.aln",
    output:
        results + "hao.gtdb.refseq.hq_anammox.remove_long_branch_1.trimal.aln",
    conda:
        "./envs/trimal.yaml"
    shell:
        "trimal -automated1 -in {input} -out {output}"

rule fasttree_after_long_branch_remove_1:
    input:
        results + "hao.gtdb.refseq.hq_anammox.remove_long_branch_1.trimal.aln",
    output:
        results + "hao.gtdb.refseq.hq_anammox.remove_long_branch_1.trimal.fasttree",
    conda:
        "./envs/fasttree.yaml"
    shell:
        "fasttree -lg {input} > {output}"

rule iqtree:
    input:
        results + "hao.gtdb.refseq.hq_anammox.remove_long_branch_1.trimal.aln"
    output:
        results + "hao.gtdb.refseq.hq_anammox.remove_long_branch_1.trimal.treefile"
    params:
        pre=results + "hao.gtdb.refseq.hq_anammox.remove_long_branch_1.trimal"
    conda:
        "./envs/iqtree.yaml"
    threads:
        12
    shell:
        "iqtree2 -s {input} -mset LG -B 1000 -nt {threads} -pre {params.pre} -redo"

"""
According to Klotz et. al. 2008 the HAO proteins have evolved from pentaheme nitrite reductase (NrfA).
Use the same sequences as they have used to root the HAO phylogeny.
In structural comparisons in the paper they find that the five hemes in the in NrfA is overlapping with
the terminal five hemes in the HAO proteins. To make this alignment we most use this haem-motifs to
anchor the alignment.
"""
rule download_root_sequences:
    input:
        accessions="../data/hao_phylogeny/outgroups.txt"
    output:
        results + "hao-outgroup.faa"
    conda:
        "./envs/entrez.yaml"
    shell:
       "epost -input {input.accessions} -db protein | efetch -format fasta > {output}"

rule merge_outgroup_hao:
    input:
        outgroup=results + "hao-outgroup.faa",
        ingroup=results + "hao.gtdb.refseq.hq_anammox.remove_long_branch_1.faa"
    output:
        results + "hao.gtdb.refseq.hq_anammox.ougtroup.faa"
    shell:
        "cat {input.ingroup} {input.outgroup} > {output}"

rule align_with_outgroup:
    input:
        results + "hao.gtdb.refseq.hq_anammox.ougtroup.faa"
    output:
        results + "hao.gtdb.refseq.hq_anammox.ougtroup.aln"
    conda:
        "./envs/mafft.yaml"
    threads:
        12
    shell:
        "mafft-linsi --reorder --thread {threads} {input} > {output}"

rule trim_with_outgroup:
    input:
        results + "hao.gtdb.refseq.hq_anammox.ougtroup.aln"
    output:
        results + "hao.gtdb.refseq.hq_anammox.ougtroup.trimal.aln"
    conda:
        "./envs/trimal.yaml"
    shell:
        "trimal -automated1 -in {input} -out {output}"

rule fasttree_with_outgroup:
    input:
        results + "hao.gtdb.refseq.hq_anammox.ougtroup.trimal.aln"
    output:
        results + "hao.gtdb.refseq.hq_anammox.ougtroup.trimal.fasttree"
    conda:
        "./envs/fasttree.yaml"
    shell:
        "fasttree -lg {input} > {output}"

rule iqtree_with_outgroup:
    input:
        results + "hao.gtdb.refseq.hq_anammox.ougtroup.trimal.aln"
    output:
        results + "hao.gtdb.refseq.hq_anammox.ougtroup.trimal.treefile"
    params:
        pre=results + "hao.gtdb.refseq.hq_anammox.ougtroup.trimal"
    threads:
        12
    conda:
        "./envs/iqtree.yaml"
    shell:
        "iqtree2 -s {input} -mset LG -nt {threads} -B 1000 -pre {params.pre} -redo"

rule align_outgroup:
    input:
        results + "hao-outgroup.faa"
    output:
        results + "hao-outgroup.aln"
    conda:
        "./envs/mafft.yaml"
    shell:
        "mafft-linsi {input} > {output}"

"""
Based on visual inspection of the alignment of the HAO-sequences and the NrfA-sequences
we can insert multiple Xs around the haem-motifs that would be used as anchors.
"""
rule anchor_ingroup:
    input:
        results+"hao.gtdb.refseq.hq_anammox.remove_long_branch_1.aln"
    output:
        results + "hao.gtdb.refseq.hq_anammox.remove_long_branch_1.anchors.faa"
    conda:
        "./envs/biopython.yaml"
    shell:
        "python add_x_ingroup.py {input} > {output}"

rule anchor_outgroup:
    input:
        results+"hao-outgroup.aln"
    output:
        results+"hao-outgroup.anchors.faa"
    conda:
        "./envs/biopython.yaml"
    shell:
        "python add_x_outgroup.py {input} > {output}"

rule merge_anchored_sequences:
    input:
        ingroup= results + "hao.gtdb.refseq.hq_anammox.remove_long_branch_1.anchors.faa",
        outgroup=results+"hao-outgroup.anchors.faa"
    output:
        results+"hao-anchors.faa"
    shell:
        "cat {input.ingroup} {input.outgroup} > {output}"

rule align_anchored_sequences:
    input:
        results+"hao-anchors.faa"
    output:
        results+"hao-anchors.aln"
    conda:
        "./envs/mafft.yaml"
    threads:
        12
    shell:
        "mafft-linsi --reorder --thread {threads} {input} > {output}"

rule remove_anchors:
    input:
        results+"hao-anchors.aln"
    output:
        results+"hao-anchors.wo_x.aln"
    conda:
        "./envs/biopython.yaml"
    shell:
        "python remove_x.py {input} > {output}"

rule trim_anchored:
    input:
        results+"hao-anchors.wo_x.aln"
    output:
        results+"hao-anchors.wo_x.trimal.aln"
    conda:
        "./envs/trimal.yaml"
    shell:
        "trimal -automated1 -in {input} -out {output}"

rule fasttree_anchored:
    input:
        results+"hao-anchors.wo_x.trimal.aln"
    output:
        results+"hao-anchors.wo_x.fasttree"
    conda:
        "./envs/fasttree.yaml"
    shell:
        "fasttree -lg {input} > {output}"

"""
Or as suggested by from 2022. The HAO/HDH can be rooted by the ONR class
of MHC proteins.
"""
rule download_onr_outgroup_query:
    output:
        results + "onr-outroup-query.faa"
    conda:
        "./envs/entrez.yaml"
    shell:
        "efetch -id AGA31970.1,ADV18468.1 -db protein -format fasta > {output}"

rule blast_onr:
    input:
        results + "onr-outroup-query.faa"
    output:
        results + "onr-outgroup-query.gtdb.blastp.tsv.gz"
    params:
        db="../data/gtdb_representatives_database/gtdb_representatives.dmnd"
    conda:
        "envs/diamond.yaml"
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
rule onr_add_accessions:
    input:
        blast=results + "onr-outgroup-query.gtdb.blastp.tsv.gz",
        prot="../data/gtdb_representatives_database/prot2accession.tsv.gz"
    output:
        results + "onr-outgroup-query.gtdb.w_accession.blastp.tsv.gz"
    shell:
        "python hzs-outside-anammox-merge-blast-accessions.py {input.blast} {input.prot} {output}"

rule onr_add_taxa:
    input:
        blast=results + "onr-outgroup-query.gtdb.w_accession.blastp.tsv.gz",
        gtdb="../data/gtdb_representatives_database/gtdb_representatives.sample_gtdb.metadata.tsv",
    output:
        results + "onr-outgroup-query.gtdb.w_taxa.blastp.tsv.gz"
    shell:
        "python hzs-outside-anammox-merge-blast-taxa.py {input.blast} {input.gtdb} {output}"

rule extract_onr_homologs:
    "Only use ADV18468.1 because the only query picks up protein containing 5 hemes"
    input:
        results + "onr-outgroup-query.gtdb.w_taxa.blastp.tsv.gz"
    output:
        results + "onr-outgroup-query.refseq.faa"
    shell:
        """
        python extract-blast-hits.py \
            --blast {input} \
            --queries ADV18468.1 \
            --scover 0.8 \
            --qcover 0.8 \
            --source refseq \
            --top 50 \
            --fasta {output}
        """

rule filter_cxxch:
    "Filter hits to contain eight hemes CxxC[H,K]"
    input:
        results + "onr-outgroup-query.refseq.faa"
    output:
        table=results + "onr-outgroup-query.refseq.cxxch_filter.tsv",
        fasta=results + "onr-outgroup-query.refseq.cxxch_filter.faa"
    conda:
        "envs/biopython.yaml"
    shell:
        "python check_cxxch.py {input} {output.table} {output.fasta}"

rule merge_onr_hao:
    input:
        hao=results + "hao.gtdb.refseq.hq_anammox.remove_long_branch_1.faa",
        onr=results + "onr-outgroup-query.refseq.cxxch_filter.faa"
    output:
        results + "hao.onr.faa"
    shell:
        "cat {input.hao} {input.onr} > {output}"

rule align_onr_hao:
    input:
        results + "hao.onr.faa"
    output:
        results + "hao.onr.aln"
    conda:
        "./envs/mafft.yaml"
    threads:
        12
    shell:
        "mafft-linsi --reorder --thread {threads} {input} > {output}"

rule trim_onr_hao:
    input:
        results + "hao.onr.aln"
    output:
        results + "hao.onr.trimal.aln"
    conda:
        "./envs/trimal.yaml"
    shell:
        "trimal -automated1 -in {input} -out {output}"

rule fasttree_onr_hao:
    input:
        results + "hao.onr.trimal.aln"
    output:
        results + "hao.onr.trimal.fasttree"
    conda:
        "envs/fasttree.yaml"
    shell:
        "fasttree -lg {input} > {output}"

rule iqtree_onr_hao:
    input:
        results + "hao.onr.trimal.aln"
    output:
        results + "hao.onr.trimal.treefile"
    params:
        pre=results + "hao.onr.trimal"
    threads:
        12
    conda:
        "./envs/iqtree.yaml"
    shell:
        "iqtree2 -s {input} -m MFP -mset WAG,LG -B 1000 -alrt 1000 -pre {params.pre} -nt {threads} -redo"

#rule
rule blast_table:
    input:
        results + "hao.gtdb.w_taxa.tsv.gz"
    output:
        results + "hao.blast.tsv"
    shell:
        "python hao_blast_table.py {input} {output}"

rule add_tree_annotation_column:
    input:
        local_taxonomy = "../data/local_taxonomy.tsv",
        local_seq = results + "hao.hq.faa",
        master=results + "hao.blast.tsv"
    output:
        results + "hao.blast.annotation.tsv"
    conda:
        "./envs/biopython.yaml"
    shell:
        "python merge_master_local_taxa.py {input.master} {input.local_taxonomy} {input.local_seq} {output}"
