backup_dir = "/media/argos-emiha442/emiha442/1_projects/1_221123_anammox_pathway_evo/processed_data/hzs-outside-anammox"
results = "../processed_data/hzs-outside-anammox"

rule all:
    input:
        #results + "/hzs.test.gtdb.tsv.gz"
        #results + "/hzs.gtdb.w_taxa.tsv.gz",
        #backup_dir + "/hzs.co-occuring.tsv",
        #backup_dir + "/hzs.co-localized.tsv",
        #results + "/hzs_a.co-localized.faa",
        #results + "/hzs_a.co-occuring.w_query.trimal.treefile",
        #results + "/hzs_a.co-occuring.w_query.trimal.fasttree",
        #results + "/hzs_bc.co-occuring.w_query.trimal.treefile",
        #results + "/hzs_a.co-occuring.w_query.trimal.aln",
        #results + "/hzs_a.co-occuring.w_query.dmnd",
        #expand(backup_dir + "/hzs_{subunit}.co-occuring.network.tsv", subunit=["a","bc"]),
        expand(results + "/{subunit}.cxxch.tsv", subunit=["hzs_a", "hzs_bc"]),
        expand(results + "/{subunit}.gtdb.hits.clustered.aln", subunit=["hzs_a", "hzs_bc"])

rule download_queries:
    output:
        results + "/query.faa"
    conda:
        "envs/entrez.yaml"
    shell:
        "efetch -id QII12198.1,QII12199.1,QII12200.1,GAX62881.1,GAX62882.1,OHC04732.1 -db protein -format fasta > {output}"

rule blast_gtdb:
    input:
        results + "/query.faa"
    output:
        results + "/hzs.gtdb.tsv.gz"
    params:
        db="../data/gtdb_representatives_database/gtdb_representatives.dmnd"
    conda:
        "envs/diamond.yaml"
    threads:
        24
    shell:
        """
        diamond blastp --db {params.db} \
            --query {input} \
            --outfmt 6 qseqid sseqid pident positive mismatch gapopen qlen slen length qstart qend sstart send evalue bitscore stitle full_sseq \
            --out {output} \
            --ultra-sensitive \
            --max-target-seqs 5000 \
            --evalue 1e-6 \
            --query-cover 50 \
            --subject-cover 50 \
            --compress 1 \
            --threads {threads}
        """

rule blast_gtdb_test:
    input:
        results + "/query.faa"
    output:
        results + "/hzs.test.gtdb.tsv.gz"
    params:
        db="../data/gtdb_representatives_database/gtdb_representatives.dmnd"
    conda:
        "envs/diamond.yaml"
    threads:
        24
    shell:
        """
        diamond blastp --db {params.db} \
            --query {input} \
            --outfmt 6 qseqid sseqid pident positive mismatch gapopen qlen slen length qstart qend sstart send evalue bitscore stitle full_sseq \
            --out {output} \
            --ultra-sensitive \
            --max-target-seqs 5000 \
            --evalue 1e-6 \
            --compress 1 \
            --threads {threads}
        """

rule add_accessions:
    input:
        blast=results + "/hzs.gtdb.tsv.gz",
        prot="../data/gtdb_representatives_database/prot2accession.tsv.gz"
    output:
        results + "/hzs.gtdb.w_accessions.tsv.gz"
    shell:
        "python hzs-outside-anammox-merge-blast-accessions.py {input.blast} {input.prot} {output}"

rule add_taxa:
    input:
        blast=results + "/hzs.gtdb.w_accessions.tsv.gz",
        gtdb="../data/gtdb_representatives_database/gtdb_representatives.sample_gtdb.metadata.tsv",
    output:
        results + "/hzs.gtdb.w_taxa.tsv.gz"
    shell:
        "python hzs-outside-anammox-merge-blast-taxa.py {input.blast} {input.gtdb} {output}"

rule hzs_co_occuring:
    input:
        results + "/hzs.gtdb.w_taxa.tsv.gz"
    output:
        co_occuring = results + "/hzs.co-occuring.tsv",
        co_localized = results + "/hzs.co-localized.tsv",
    shell:
        "python hzs-outside-anammox-find-colocalized.py {input} {output.co_occuring} {output.co_localized}"

rule hzs_a_queries:
    input:
        query = results + "/query.faa"
    output:
        results + "/query_hzs_a.faa"
    conda:
        "envs/seqtk.yaml"
    shell:
        """
        echo "GAX62882.1,QII12200.1" | tr "," "\\n" | seqtk subseq {input} - > {output}
        """

rule hzs_bc_queries:
    input:
        query = results + "/query.faa"
    output:
        results + "/query_hzs_bc.faa"
    conda:
        "envs/seqtk.yaml"
    shell:
        """
        echo "GAX62881.1,QII12198.1,QII12199.1" | tr "," "\\n" | seqtk subseq {input} - > {output}
        """
rule extract_all_hzs_hits:
    input:
        results + "/hzs.gtdb.w_taxa.tsv.gz"
    output:
        results + "/{subunit}.gtdb.hits.faa"
    params:
        subunit = "{subunit}"
    shell:
        "python hzs-outside-anammox-extract-all.py {input} {params.subunit} > {output}"

rule cluster_hits:
    input:
        results + "/{subunit}.gtdb.hits.faa"
    output:
        results + "/{subunit}.gtdb.hits.clustered.faa"
    conda:
        "envs/cd-hit.yaml"
    shell:
       "cd-hit -i {input} -o {output} -c 0.8"

rule align_all_hzs_hits:
    input:
        results + "/{subunit}.gtdb.hits.clustered.faa"
    output:
        results + "/{subunit}.gtdb.hits.clustered.aln"
    conda:
        "./envs/mafft.yaml"
    threads:
        12
    shell:
        "mafft-linsi --thread {threads} {input} > {output}"        
    

rule get_hit_accessions:
    input:
        results + "/hzs.gtdb.w_taxa.tsv.gz"
    output:
        results + "/{subunit}.gtdb.hits.accessions"
    params: 
        subunit = "{subunit}"
    shell:
        "python hzs-outside-anammox-extract-accessions.py {input} {params.subunit} > {output}" 

rule get_hit_genome_info:
    input:
        results + "/{subunit}.gtdb.hits.accessions"
    output:
        results + "/{subunit}.gtdb.hits.genome_info.tsv"
    conda:
        "./envs/ncbi-datasets.yaml"
    shell:
        """
        datasets summary genome accession --inputfile {input} --as-json-lines | \
        dataformat tsv genome > {output}
        """

rule merge_blast_and_structure:
    input:
        blast=results+"/hzs.gtdb.w_taxa.tsv.gz",
        structure_alignment="../processed_data/hzs-structure-alignment/alignment_summary.tsv",
        structure_map="../processed_data/hzs-structure-alignment/{subunit}_alphafold_identifiers.tsv"
    output:
        results+"/{subunit}.blast.structure.tsv"
    params:
        subunit="{subunit}"
    shell:
       "python merge_hzs_blast_structur_data.py {input.blast} {input.structure_map} {input.structure_alignment} {output} {params.subunit}"

rule add_MAG_info:
    input:
        table="/{subunit}.blast.structure.tsv"
        genome_info=results + "/{subunit}.gtdb.hits.genome_info.tsv"
    output:
        results + "/{subunit}.blast.structure.genome.tsv"
    shell:
        "python add_MAG_info.py {input.table} {input.genome_info} {output}

    

rule extract_co_occuring_hzs_a:
    """
    Extract hits to HZS-A from species that also have a hit to HZS-BC
    """
    input:
        results + "/hzs.co-occuring.tsv"
    output:
        results + "/hzs_a.co-occuring.faa"
    shell:
        "python hzs-outside-anammox-extract-hzs-a.py {input} > {output}"

rule extract_co_occuring_hzs_bc:
    """
    Extract hits to HZS-BC from species that also have a hit to HZS-BC
    """
    input:
        results + "/hzs.co-occuring.tsv"
    output:
        results + "/hzs_bc.co-occuring.faa"
    shell:
        "python hzs-outside-anammox-extract-hzs-bc.py {input} > {output}"

rule merge_co_occuring_hzs_a_w_query:
    input:
        hzs_a = results + "/hzs_a.co-occuring.faa",
        query = results + "/query_hzs_a.faa"
    output:
        results + "/hzs_a.co-occuring.w_query.faa"
    shell:
        "cat {input.hzs_a} {input.query} > {output}"

rule merge_co_occuring_hzs_bc_w_query:
    input:
        hzs_bc = results + "/hzs_bc.co-occuring.faa",
        query = results + "/query_hzs_bc.faa"
    output:
        results + "/hzs_bc.co-occuring.w_query.faa"
    shell:
        "cat {input.hzs_bc} {input.query} > {output}"

rule align_co_occuring_hzs_a:
    input:
        results + "/hzs_a.co-occuring.w_query.faa"
    output:
        results + "/hzs_a.co-occuring.w_query.aln"
    conda:
        "envs/mafft.yaml"
    threads:
        12
    shell:
        "mafft-linsi --reorder --thread {threads} {input} > {output}"

rule align_co_occuring_hzs_bc:
    input:
        results + "/hzs_bc.co-occuring.w_query.faa"
    output:
        results + "/hzs_bc.co-occuring.w_query.aln"
    conda:
        "envs/mafft.yaml"
    threads:
        12
    shell:
        "mafft-linsi --reorder --thread {threads} {input} > {output}"

rule trim_co_occuring_hzs_a_alignment:
    input:
        results + "/hzs_a.co-occuring.w_query.aln"
    output:
        results + "/hzs_a.co-occuring.w_query.trimal.aln"
    conda:
        "envs/trimal.yaml"
    shell:
        "trimal -gt 0.5 -in {input} -out {output}"

rule trim_co_occuring_hzs_bc_alignment:
    input:
        results + "/hzs_bc.co-occuring.w_query.aln"
    output:
        results + "/hzs_bc.co-occuring.w_query.trimal.aln"
    conda:
        "envs/trimal.yaml"
    shell:
        "trimal -gt 0.5 -in {input} -out {output}"

rule fasttree_co_occuring_hzs_a:
    input:
        results + "/hzs_a.co-occuring.w_query.trimal.aln"
    output:
        results + "/hzs_a.co-occuring.w_query.trimal.fasttree"
    conda:
        "envs/fasttree.yaml"
    shell:
        "fasttree -lg {input} > {output}"

rule iqtree_co_occuring_hzs_a:
    input:
        results + "/hzs_a.co-occuring.w_query.trimal.aln"
    output:
        results + "/hzs_a.co-occuring.w_query.trimal.treefile"
    params:
        prefix=results + "/hzs_a.co-occuring.w_query.trimal"
    conda:
        "envs/iqtree.yaml"
    threads:
        12
    shell:
        "iqtree2 -s {input} -m LG+F+G4 -bb 1000 -nt {threads} -pre {params.prefix}"

rule iqtree_co_occuring_hzs_bc:
    input:
        results + "/hzs_bc.co-occuring.w_query.trimal.aln"
    output:
        results + "/hzs_bc.co-occuring.w_query.trimal.treefile"
    params:
        prefix=results + "/hzs_bc.co-occuring.w_query.trimal"
    conda:
        "envs/iqtree.yaml"
    threads:
        12
    shell:
        "iqtree2 -s {input} -m LG+F+G4 -bb 1000 -nt {threads} -pre {params.prefix} -redo"


rule extract_co_localized_hzs_a:
    input:
        results + "/hzs.co-localized.tsv",
    output:
        results + "/hzs_a.co-localized.faa"
    shell:
        "python ./hzs-outside-anammox-extract-hzs-a.py {input} > {output}"


rule merge_hzs_a_w_query:
    input:
        hzs_a = results + "/hzs_a.co-localized.faa",
        query = results + "/query_hzs_a.faa"
    output:
        results + "/hzs_a.co-localized.w_query.faa"
    shell:
        "cat {input.hzs_a} {input.query} > {output}"

# Create a protein similarity network by doing all against all blast
rule hzs_co_occuring_db:
    input:
        results + "/hzs_{subunit}.co-occuring.w_query.faa"
    output:
        results + "/hzs_{subunit}.co-occuring.w_query.dmnd"
    params:
        prefix = results + "/hzs_{subunit}.co-occuring.w_query"
    conda:
        "envs/diamond.yaml"
    shell:
        "diamond makedb --db {params.prefix} --in {input}"


rule hzs_subunit_all_vs_all:
    input:
        query=results + "/hzs_{subunit}.co-occuring.w_query.faa",
        db=results + "/hzs_{subunit}.co-occuring.w_query.dmnd"
    output:
        results + "/hzs_{subunit}.co-occuring.w_query.all_vs_all.blastp.tsv"
    conda:
        "envs/diamond.yaml"
    threads:
        24
    shell:
        """
        diamond blastp --db {input.db} \
            --query {input.query} \
            --outfmt 6 qseqid sseqid pident ppos positive mismatch gapopen qlen slen length qstart qend sstart send evalue bitscore \
            --out {output} \
            --ultra-sensitive \
            --max-target-seqs 5000 \
            --evalue 1e-6 \
            --query-cover 50 \
            --subject-cover 50 \
            --threads {threads}
        """

rule hzs_subunit_protein_similiarity_network:
    """
    Output of this rule is imported to cytoscape for visualization
    """
    input:
        original_blast = "../processed_data/hzs-outside-anammox/hzs.gtdb.w_taxa.tsv.gz",
        all_vs_all_blast = "../processed_data/hzs-outside-anammox/hzs_{subunit}.co-occuring.w_query.all_vs_all.blastp.tsv",
        co_localized = "../processed_data/hzs-outside-anammox/hzs.co-localized.tsv"
    output:
        "../processed_data/hzs-outside-anammox/hzs_{subunit}.co-occuring.network.tsv"
    shell:
        "python hzs-outside-anammox-create-network-data.py {input.original_blast} {input.all_vs_all_blast} {input.co_localized} {output}"

rule search_cxxch_motifs:
    input:
        "../processed_data/hzs-outside-anammox/{subunit}.gtdb.hits.faa"
    output:
        "../processed_data/hzs-outside-anammox/{subunit}.cxxch.tsv"
    conda:
        "./envs/biopython.yaml"
    shell:
        "python hzs-outside-anammox-cxxch.py {input} {output}"

# Copy files to Argos
rule backup_hzs_co_occuring:
    input:
        results + "/hzs.co-occuring.tsv"
    output:
        backup_dir + "/hzs.co-occuring.tsv"
    shell:
        "cp {input} {output}"

rule backup_hzs_co_localized:
    input:
        results + "/hzs.co-localized.tsv"
    output:
        backup_dir + "/hzs.co-localized.tsv"
    shell:
        "cp {input} {output}"
rule backup_network:
    input:
        "../processed_data/hzs-outside-anammox/hzs_{subunit}.co-occuring.network.tsv"
    output:
        backup_dir + "/hzs_{subunit}.co-occuring.network.tsv"
    shell:
        "cp {input} {output}"
