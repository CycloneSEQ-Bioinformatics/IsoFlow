rule count_reads:
    input:
        bam = rules.map_reads_to_transciptome.output.sbam,
        trs = rules.build_transcriptome_index.output.trs_fa
    output:
        tsv = "count_reads/{sample}/quant.sf",
    params:
        tsv_dir = "count_reads/{sample}",
        libtype = config["salmon_libtype"]
    log:
        "logs/count_reads.{sample}.log"
    benchmark:
        "benchmarks/count_reads.{sample}.benchmark"
    threads: config["threads"]["count_reads"]
    shell: """
    salmon quant --noErrorModel -p {threads} -t {input.trs} -l {params.libtype} -a {input.bam} -o {params.tsv_dir} 2>{log}
    """

rule merge_counts:
    input:
        count_tsvs = expand("count_reads/{sample}/quant.sf", sample=all_samples.keys()),
        trs_info = rules.build_transcriptome_index.output.trs_info
    output:
        count_tsv = "merge_counts/all_counts.tsv",
        tpm_tsv = "merge_counts/all_TPM.tsv",
    benchmark:
        "benchmarks/merge_counts.benchmark"
    shell:"""
    {SNAKEDIR}/scripts/merge_count_tsvs.py -z -f transcript_id -o {output.count_tsv} -tsvs {input.count_tsvs}
    {SNAKEDIR}/scripts/merge_count_tsvs.py -z -f transcript_id -tpm True -o {output.tpm_tsv} -tsvs {input.count_tsvs}
    {SNAKEDIR}/scripts/add_gene_name.R {input.trs_info} {output.count_tsv} transcript_id
    {SNAKEDIR}/scripts/add_gene_name.R {input.trs_info} {output.tpm_tsv} transcript_id
    """
