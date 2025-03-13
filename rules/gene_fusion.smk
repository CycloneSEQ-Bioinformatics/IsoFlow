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

rule plot_fusions:
    input:
        csvs = expand("jaffal_fusion/{sample}/jaffa_results.csv", sample=all_samples.keys())
    output:
        cat_csv = "jaffal_fusion/jaffa_results.csv",
        plot_log = "jaffal_fusion/plot.log",
    params:
        ref_version = config["reference"][config["specie"]]["version"],
        min_span_reads = 30,
        outdir = 'jaffal_fusion',
    log: "logs/plot_fusions.log"
    benchmark: "benchmarks/plot_fusions.log"
    shell:"""
    awk 'FNR==1 && NR!=1{{next}} {{print}}' {input.csvs} | sed 's/.usable.fq//' > {params.outdir}/jaffa_results.csv
    {SNAKEDIR}/scripts/plot_fusions.R {params.outdir}/jaffa_results.csv {params.min_span_reads} {params.ref_version} {params.outdir} > {log} 2>&1
    echo "Plot successfully" > {output.plot_log}
    """
