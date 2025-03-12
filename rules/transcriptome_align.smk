rule build_transcriptome_index:
    input:
        trs_gtf = config["reference"][config["specie"]]["gtf"],
        genome = config["reference"][config["specie"]]["fasta"]
    output:
        trs_fa = "transcriptome_index/transcripts.fa",
        index = "transcriptome_index/transcriptome_index.mmi",
        trs_info = "transcriptome_index/transcripts_info.tsv"
    params:
        with_ERCC = config["with_ERCC"],
        ERCC_fasta = config["reference"]["ERCC"]["fasta"],
    log:
        "logs/transcriptome_index.log"
    benchmark:
        "benchmarks/transcriptome_index.benchmark"
    threads: config["threads"]["build_transcriptome_index"]
    shell:"""
    gffread -w {output.trs_fa} -g {input.genome} {input.trs_gtf}
    if [ {params.with_ERCC} == "True" ];then
        cat {params.ERCC_fasta} >> {output.trs_fa}
    fi
    minimap2 -t {threads} -d {output.index} {output.trs_fa} 2>{log}
    {SNAKEDIR}/scripts/make_trx_infoTab.R {input.trs_gtf} {output.trs_info} 2>>{log}
    """

rule map_reads_to_transciptome:
    input:
        index = rules.build_transcriptome_index.output.index,
        fastq = "fl_reads/{sample}/{sample}.usable.fq.gz"
    output:
        bam = "transcriptome_alignments/{sample}/aligned.bam",
        sbam = "transcriptome_alignments/{sample}/aligned.sort.bam",
        alnstat = "transcriptome_alignments/{sample}/read_aln_stats.tsv"
    params:
        opts = config["minimap2_opts"]["transcriptome"]
    log:
        "logs/map_reads_to_transciptome.{sample}.log"
    benchmark:
        "benchmarks/map_reads_to_transciptome.{sample}.benchmark"
    threads: config["threads"]["map_reads"]
    shell:"""
    minimap2 -t {threads} {params.opts} {input.index} {input.fastq} 2>{log} \
    | samtools view -Sb > {output.bam}
    samtools sort -@ {threads} --write-index -o {output.sbam} {output.bam} 2>>{log}
    samtools flagstat -@ {threads} -O tsv {output.sbam} | awk -f {SNAKEDIR}/scripts/bamstat_fmt.awk > {output.alnstat}
    """

rule calc_gene_trx_cov:
    input:
        bam = "transcriptome_alignments/{sample}/aligned.sort.bam",
        trs_info = rules.build_transcriptome_index.output.trs_info
    output:
        cov = "known_transcripts_coverage/{sample}/coverage.tsv"
    params:
        outdir =  "known_transcripts_coverage/{sample}/"
    log: "logs/calc_gene_trx_cov.{sample}.log"
    benchmark: "benchmarks/calc_gene_trx_cov.{sample}.benchmark"
    threads: 4
    shell:"""
    mosdepth -n -t {threads} {params.outdir}/out {input.bam} > {log} 2>&1
    {SNAKEDIR}/scripts/get_gene_trx_cov.R {input.trs_info} {params.outdir}/out.mosdepth.global.dist.txt {output.cov} >> {log} 2>&1
    """

rule plot_gene_trx_cov:
    input:
        covs = expand("known_transcripts_coverage/{sample}/coverage.tsv", sample=all_samples.keys())
    output:
        cov = "known_transcripts_coverage/coverage.tsv",
        png = "known_transcripts_coverage/coverage.png"
    params:
        preffix = "known_transcripts_coverage/coverage"
    log: "logs/plot_gene_trx_cov.log"
    benchmark: "benchmarks/plot_gene_trx_cov.benchmark"
    shell:"""
    {SNAKEDIR}/scripts/plot_gene_trx_cov.R {input.covs} {params.preffix} > {log} 2>&1
    """
