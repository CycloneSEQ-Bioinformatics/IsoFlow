rule isoquant:
    input:
        bams = expand("genome_alignments/{sample}/aligned.sort.bam", sample=all_samples.keys()),
        genome = config["reference"][config["specie"]]["fasta"],
        gtf = config["reference"][config["specie"]]["gtf"]
    output:
        trs_gtf = "isoquant/OUT/OUT.transcript_models.gtf",
        gene_counts = "isoquant/OUT/OUT.gene_counts.tsv",
        transcript_counts = "isoquant/OUT/OUT.transcript_counts.tsv"
    params:
        samples = expand("{sample}", sample=all_samples.keys()),
        outdir = "isoquant",
        stranded = config["stranded"]
    log:
        "logs/isoquant.log"
    benchmark:
        "benchmarks/isoquant.benchmark"
    threads: config["threads"]["isoquant"]
    shell:"""
    isoquant.py --reference {input.genome} --genedb {input.gtf} --complete_genedb --bam {input.bams} --labels {params.samples} --data_type nanopore --stranded {params.stranded} -o {params.outdir} --threads {threads} 2>{log}
    """

#rule run_gffcompare:
#    input:
#        trs_gtf = rules.isoquant.output.trs_gtf,
#        annotation_gtf = config["reference"][config["specie"]]["gtf"],
#    output:
#        tracking = "gffcompare/transcript_models.tracking"
#    params:
#        outdir = "gffcompare",
#        preffix = "transcript_models"
#    log: "logs/run_gffcompare.log"
#    benchmark: "benchmarks/run_gffcompare.benchmark"
#    shell:"""
#    gffcompare -o {params.outdir}/{params.preffix} -r {input.annotation_gtf} -R {input.trs_gtf} > {log} 2>&1
#    """

rule sqanti3_qc:
    input:
        trs_gtf = rules.isoquant.output.trs_gtf,
        annotation_gtf = config["reference"][config["specie"]]["gtf"],
        genome = config["reference"][config["specie"]]["fasta"]
    output:
        pdf = "sqanti_qc/OUT.transcript_models_SQANTI3_report.pdf",
        classification = "sqanti_qc/OUT.transcript_models_classification.txt",
        junction = "sqanti_qc/OUT.transcript_models_junctions.txt",
        png = "sqanti_qc/OUT.transcript_models.evidence.png"
    params:
        SQANTI3_HOME = config["SQANTI3_HOME"],
        outdir = "sqanti_qc"
    log:
        "logs/sqanti_qc.log"
    benchmark:
        "benchmarks/sqanti_qc.benchmark"
    threads: config["threads"]["sqanti3_qc"]
    conda: "SQANTI3.env"
    shell:"""
    python {params.SQANTI3_HOME}/sqanti3_qc.py {input.trs_gtf} {input.annotation_gtf} {input.genome} --CAGE_peak {params.SQANTI3_HOME}/data/ref_TSS_annotation/human.refTSS_v3.1.hg38.bed --polyA_motif_list {params.SQANTI3_HOME}/data/polyA_motifs/mouse_and_human.polyA_motif.txt -d {params.outdir} --cpus {threads} --report pdf > {log} 2>&1
    {SNAKEDIR}/scripts/plot_sqanti3_qc.R {output.classification} {output.png} > {log} 2>&1
    """
