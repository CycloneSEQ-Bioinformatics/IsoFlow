rule full_len_reads:
    input:
        pass_fastq = lambda wildcards: all_samples[wildcards.sample]
    output:
        fullLen_fastq = "fl_reads/{sample}/{sample}.full-length.fq.gz",
        rescue_fastq = "fl_reads/{sample}/{sample}.rescued.fq.gz",
        merge_fastq = "fl_reads/{sample}/{sample}.usable.fq.gz"
    params:
        sample_name = "{sample}",
        opts = config["glycine_opts"],
        outdir = "fl_reads/{sample}"
    log:
        "logs/full_len_reads.{sample}.log"
    benchmark:
        "benchmarks/full_len_reads.{sample}.benchmark"
    threads: config["threads"]["glycine"]
    shell:"""
    glycine -f {input.pass_fastq} -n {params.sample_name} -t {threads} {params.opts} -o {params.outdir} 2>{log}
    cat {output.fullLen_fastq} {output.rescue_fastq} > {output.merge_fastq}
    """

rule qc_pass:
    input:
        fastq = lambda wildcards: all_samples[wildcards.sample],
    output:
        stat = "fl_reads/{sample}/{sample}.passed.stat.gz",
        plot_read_len = "fl_reads/{sample}/{sample}.passed.length_distribution.png",
        plot_read_gc = "fl_reads/{sample}/{sample}.passed.gc_content.png",
        plot_read_qual = "fl_reads/{sample}/{sample}.passed.quality_distribution.png"
    params:
        preffix = "fl_reads/{sample}/{sample}.passed",
    log:
        "logs/qc_pass.{sample}.log"
    benchmark:
        "benchmarks/qc_pass.{sample}.benchmark"
    threads: 8
    shell:"""
    if [ ! -f {params.preffix}.sampled.fq.gz ];then
        seqkit sample -p 0.1 -j {threads} -o {params.preffix}.sampled.fq.gz {input.fastq} > {log} 2>&1
    fi
    if [ ! -f {output.stat} ];then
        seqkit fx2tab -l -g -q -n -j {threads} --compress-level 5 -o {output.stat} {params.preffix}.sampled.fq.gz > {log} 2>&1
    fi
    {SNAKEDIR}/scripts/plot_readQC.py {output.stat} "passed reads" {params.preffix} > {log} 2>&1
    """

rule qc_fulllen:
    input:
        fastq = "fl_reads/{sample}/{sample}.usable.fq.gz"
    output:
        stat = "fl_reads/{sample}/{sample}.usable.stat.gz",
        plot_read_len = "fl_reads/{sample}/{sample}.usable.length_distribution.png",
        plot_read_gc = "fl_reads/{sample}/{sample}.usable.gc_content.png",
        plot_read_qual = "fl_reads/{sample}/{sample}.usable.quality_distribution.png"
    params:
        preffix = "fl_reads/{sample}/{sample}.usable"
    log:
        "logs/qc_fullLen.{sample}.log"
    benchmark:
        "benchmarks/qc_fullLen.{sample}.benchmark"
    threads: 8
    shell:"""
    if [ ! -f {params.preffix}.sampled.fq.gz ];then
        seqkit sample -p 0.1 -j {threads} -o {params.preffix}.sampled.fq.gz {input.fastq} > {log} 2>&1
    fi
    if [ ! -f {output.stat} ];then
        seqkit fx2tab -l -g -q -n -j {threads} --compress-level 5 -o {output.stat} {params.preffix}.sampled.fq.gz > {log} 2>&1
    fi
    {SNAKEDIR}/scripts/plot_readQC.py {output.stat} "full-length reads" {params.preffix} > {log} 2>&1
    """
