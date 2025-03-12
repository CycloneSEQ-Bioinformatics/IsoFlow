rule bam_to_bigwig:
    input:
        bam = "genome_alignments/{sample}/aligned.sort.bam"
    output:
        bigwig = "known_transcripts_depth/{sample}.bw"
    log:
        "logs/bam_to_bigwig.{sample}.log"
    benchmark:
        "benchmarks/bam_to_bigwig.{sample}.benchmark"
    threads: 8
    shell:"""
    bamCoverage -b {input.bam} -o {output.bigwig} --normalizeUsing CPM -p {threads} > {log}
    """

rule get_bed12:
   input:
       gtf = config["reference"][config["specie"]]["gtf"]
   output:
       bed12 = "known_transcripts_depth/transcripts.bed12"
   shell:"""
   awk -f {SNAKEDIR}/scripts/gtf2bed12.awk {input.gtf} > {output.bed12}
   """

rule compute_matrix:
    input:
        bigwigs = expand("known_transcripts_depth/{sample}.bw", sample=all_samples.keys()),
        bed12 = rules.get_bed12.output.bed12
    output:
        readProfileMatrix = "known_transcripts_depth/readProfileMatrix.tsv.gz"
    log:
        "logs/compute_matrix.log"
    benchmark:
        "benchmarks/compute_matrix.benchmark"
    threads: 16
    shell:"""
    computeMatrix scale-regions -S {input.bigwigs} -R {input.bed12} -o {output.readProfileMatrix} --upstream 1000 --downstream 1000 --sortRegions ascend --missingDataAsZero --skipZeros --metagene -p {threads}
    """

rule plot_profile:
    input:
        readProfileMatrix = rules.compute_matrix.output.readProfileMatrix
    output:
        png = "known_transcripts_depth/Profile.png"
    log: "logs/plot_profile.log"
    benchmark: "benchmarks/plot_profile.benchmark"
    shell:"""
    plotProfile -m {input.readProfileMatrix} -o {output.png} --perGroup --plotType se --yAxisLabel 'mean CPM'
    """
