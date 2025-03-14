rule build_genome_index: ## build minimap2 genome index
    input:
        genome = config["reference"][config["specie"]]["fasta"]
    output:
        index = "genome_index/genome_splice_index.mmi"
    params:
        opts = config["minimap_genome_index_opts"]
    log: "logs/build_genome_index.log"
    benchmark: "benchmarks/build_genome_index.benchmark"
    threads: config["threads"]["build_genome_index"]
    shell:"""
    minimap2 -t {threads} {params.opts} -d {output.index} {input.genome} 2>{log}
    """

rule map_reads_to_genome: ## map reads using minimap2
    input:
       index = rules.build_genome_index.output.index,
       fastq = "fl_reads/{sample}/{sample}.usable.fq.gz"
    output:
       bam = "genome_alignments/{sample}/aligned.bam",
       sbam = "genome_alignments/{sample}/aligned.sort.bam",
       alnstat = "genome_alignments/{sample}/read_aln_stats.tsv"
    params:
        opts = config["minimap2_opts"]["genome"]
    log:
        "logs/map_reads_to_genome.{sample}.log"
    benchmark:
        "benchmarks/map_reads_to_genome.{sample}.benchmark"
    threads: config["threads"]["map_reads"]
    shell:"""
    minimap2 -t {threads} {params.opts} {input.index} {input.fastq} 2>{log} | samtools view -Sb > {output.bam}
    samtools view -h -q 10 -F 2304 {output.bam} | samtools sort -@ {threads} --write-index -o {output.sbam} 2>>{log}
    samtools flagstat -@ {threads} -O tsv {output.bam} | awk -f {SNAKEDIR}/scripts/bamstat_fmt.awk > {output.alnstat}
    """
