rule jaffal_fuison:
    input:
        fastq = "fl_reads/{sample}/{sample}.usable.fq.gz"
    output:
        csv = 'jaffal_fusion/{sample}/jaffa_results.csv',
        fasta = 'jaffal_fusion/{sample}/jaffa_results.fasta'
    params:
        refBase = config["jaffal_refBase"],
        outdir = 'jaffal_fusion/{sample}/'
    log:
        "logs/jaffal_fuison.{sample}.log"
    benchmark:
        "benchmarks/jaffal_fuison.{sample}.benchmark"
    threads: config["threads"]["jaffal"]
    conda: "jaffa_env"
    shell:"""
    groovy=$(which bpipe | xargs -n 1 dirname)/../share/jaffa-2.3-0/JAFFAL.groovy
    bpipe run -n {threads} -p refBase={params.refBase} -p jaffa_output={params.outdir} $groovy {input.fastq} >{log} 2>&1
    """
