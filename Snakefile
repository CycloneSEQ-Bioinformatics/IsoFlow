import os
import pandas as pd
from collections import OrderedDict

configfile: "config.yaml"
SNAKEDIR = os.path.dirname(workflow.snakefile)
SUPPA_HOME = config["SUPPA_HOME"]

def load_manifest(infile):
    manifest_df = pd.read_table(infile, header=None, names=["sample", "alias", "condition", "fastq"])
    manifest_dict = {}
    all_samples_dict = {}
    for row_tuple in manifest_df.itertuples():
        all_samples_dict[row_tuple.alias] = row_tuple.fastq
    return manifest_df, all_samples_dict, 

manifest, all_samples = load_manifest(config["manifest"])

# preprocess and QC
include: f"{SNAKEDIR}/rules/preprocess_qc.smk"
# module1: transcriptome reconstruction and noval isoform discovery
if config["modules"]["module1"]:
    include: f"{SNAKEDIR}/rules/genome_align.smk"
    include: f"{SNAKEDIR}/rules/known_transcripts_depth.smk"
    include: f"{SNAKEDIR}/rules/transcriptome_reconstruction.smk"
# module2: transcrips quantification and differential expression/splicing analysis
if config["modules"]["module2"]:
    include: f"{SNAKEDIR}/rules/transcriptome_align.smk"
    include: f"{SNAKEDIR}/rules/transcriptome_quant.smk"
    include: f"{SNAKEDIR}/rules/diff_exp.smk"
    include: f"{SNAKEDIR}/rules/diff_splice.smk"
# module3: gene fusion detection
if config["modules"]["module3"]:
    include: f"{SNAKEDIR}/rules/gene_fusion.smk"


def rule_all_input(modules_dict):
    all_input = list()
    all_input.extend(expand("fl_reads/{sample}/{sample}.passed.stat.gz", sample=all_samples.keys()))
    all_input.extend(expand("fl_reads/{sample}/{sample}.passed.length_distribution.png", sample=all_samples.keys()))
    all_input.extend(expand("fl_reads/{sample}/{sample}.passed.gc_content.png", sample=all_samples.keys()))
    all_input.extend(expand("fl_reads/{sample}/{sample}.passed.quality_distribution.png", sample=all_samples.keys()))
    all_input.extend(expand("fl_reads/{sample}/{sample}.usable.stat.gz", sample=all_samples.keys()))
    all_input.extend(expand("fl_reads/{sample}/{sample}.usable.length_distribution.png", sample=all_samples.keys()))
    all_input.extend(expand("fl_reads/{sample}/{sample}.usable.gc_content.png", sample=all_samples.keys()))
    all_input.extend(expand("fl_reads/{sample}/{sample}.usable.quality_distribution.png", sample=all_samples.keys()))
    if modules_dict['module1']:
        all_input.append("known_transcripts_depth/Profile.png")
        all_input.append("sqanti_qc/OUT.transcript_models_SQANTI3_report.pdf")
        all_input.append("sqanti_qc/OUT.transcript_models_classification.txt")
        all_input.append("sqanti_qc/OUT.transcript_models_junctions.txt")
        all_input.append("sqanti_qc/OUT.transcript_models.evidence.png")
        #all_input.append("gffcompare/transcript_models.tracking")
    if modules_dict['module2']:
        all_input.append("known_transcripts_coverage/coverage.png")
        all_input.append("diff_exp/results_dge.tsv")
        all_input.append("diff_exp/results_dge.maplot.png")
        all_input.append("diff_exp/results_dge.volcanoplot.png")
        all_input.append("diff_exp/results_dte.tsv")
        all_input.append("diff_exp/results_dte.maplot.png")
        all_input.append("diff_exp/results_dte.volcanoplot.png")
        all_input.append("diff_splice/plot_sashimi/plot.log")
        #all_input.append("de_analysis/results_dge.tsv")
        #all_input.append("de_analysis/results_dge.pdf")
        #all_input.append("de_analysis/results_dtu.pdf")
        #all_input.append("de_analysis/results_dtu_dexseq_gene.tsv")
        #all_input.append("de_analysis/results_dtu_dexseq_transcript.tsv")
        #all_input.append("de_analysis/results_dtu_dexseq_transcript.detail.tsv")
        #all_input.append("de_analysis/results_dtu_stageR.tsv")
        #all_input.append("de_analysis/filtered_transcript_counts_with_genes.tsv")
        #all_input.append("de_analysis/filtered_transcript_counts_with_genes.tsv")
        #all_input.append("de_analysis/dtu_plots.pdf")
    if modules_dict['module3']:
        all_input.append("jaffal_fusion/jaffa_results.csv")
        all_input.append("jaffal_fusion/plot.log")
       
    return all_input


rule all:
    input:
        "Report.html"

rule generate_report:
    input:
        rule_all_input(config["modules"]),
    output:
        html = "Report.html",
    shell:"""
    {SNAKEDIR}/scripts/ReportGeneration.py
    """

