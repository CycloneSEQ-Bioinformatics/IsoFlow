# IsoFlow: Bulk Full-Length Transcriptome Analysis Pipeline

[中文](README.CHN.md)

## Introduction

![](pictures/IsoFlow_overview.png)

This pipeline is designed for CycloneSEQ bulk full-length transcriptome analysis and is divided into three independent modules:

**Module 1**: Transcriptome reconstruction and novel isoform discovery.

**Module 2**: Transcript quantification and differential expression/splicing analysis.

**Module 3**: Gene fusion detection.

By default, the pipeline executes all three modules.

If users wish to run only specific modules, they can modify the configuration file `config.yaml` to enable or disable the corresponding modules.

## Install and use

### Install

(1) Install miniconda by following the [official guide](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html).

(2) Copy the folder of this pipeline to your local machine.

(3) Create a conda environment specifically for this analysis pipeline.
    
```shell
conda env create -n pipeline_fulllong_transcriptome -f env.yaml
```

(4) Install SUPPA2

```shell
wget -O SUPPA-2.4.tar.gz https://github.com/comprna/SUPPA/archive/refs/tags/v2.4.tar.gz
tar zxvf SUPPA-2.4.tar.gz
````

Modify the config.yaml file to provide the installation path of SUPPA.

```shell
SUPPA_HOME: "/PATH/TO/SUPPA-2.4"
```

(5) Install JAFFA

JAFFA is a tool designed for gene fusion detection. Its runtime environment is relatively unique, requiring the creation of a dedicated conda environment and the installation of the tool within this environment. The specific steps are as follows:

```shell
conda create -n jaffa_env python=3.8
conda activate jaffa_env
conda install -c bioconda bpipe=0.9.11 jaffa=2.3 minimap2=2.28

target_dir=$(which bpipe | xargs -n 1 dirname)
cp $target_dir/../share/jaffa-2.3-0/docker/tools.groovy $target_dir/../share/jaffa-2.3-0/
echo "make_3_gene_fusion_table='make_3_gene_fusion_table'" >> $target_dir/../share/jaffa-2.3-0/tools.groovy
for i in make_3_gene_fusion_table.c++ extract_seq_from_fasta.c++ make_simple_read_table.c++ process_transcriptome_align_table.c++;do g++ -std=c++11 -O3 -o $target_dir/$(basename $i .c++) $target_dir/../share/jaffa-2.3-0/src/$i;done

g++ -std=c++11 -O3 -o $target_dir/make_count_table $target_dir/../share/jaffa-2.3-0/src/make_count_table.c++
```

(6) Install SQANTI3

SQANTI3 is a tool designed for the in-depth characterization of isoforms obtained by full-length transcript sequencing. Its runtime environment is relatively unique, requiring the creation of a dedicated conda environment and the installation of the tool within this environment. The specific steps are as follows:

```bash
wget -O SQANTI3-5.2.2.tar.gz https://github.com/ConesaLab/SQANTI3/archive/refs/tags/v5.2.2.tar.gz
tar zxvf SQANTI3-5.2.2.tar.gz
cd SQANTI3-5.2.2

conda env create -f SQANTI3.conda_env.yml

# Check and test whether the conda environment has been successfully created
conda env list
conda activate SQANTI3.env
```

(7) Install ggsashimi

ggsashimi is a tool designed for visualizing differential alternative splicing events. Its runtime environment is relatively unique, requiring the creation of a dedicated conda environment and the installation of the tool within this environment. The specific steps are as follows:

```shell
conda config --set channel_priority flexible
conda env create -n ggsashimi_env -f env.ggsashimi.yaml
```

(8) Install Glycine

Glycine is a proprietary tool developed by CycloneSEQ for identifying full-length reads. It can recognize full-length reads from raw sequencing data based on the sequence composition and relative structure of dual-end amplification primers (TSO & TRP) and polyA/T.

Since the provided software is an executable file, simply save it to a directory and add it to the environment variables.

```shell
wget -c https://github.com/CycloneSEQ-Bioinformatics/Glycine/releases/download/v1.0.0/glycine.tar.gz
tar zxvf glycine.tar.gz
```   
    
### Required input files

Common input files, such as the reference genome, and configuration parameters are specified in the config.yaml file.

Preparation and Specification of Input Files

- Reference genome sequences and annotations for (multiple) target species.

    ```yaml
    reference:
        Human:
            fasta: "/data/resources/reference/hg38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
            gtf: "/data/resources/reference/hg38/Homo_sapiens.GRCh38.113.primary_assembly.gtf"
            version: "hg38"
        Mouse:
            fasta: "/data/resources/reference/GRCm39/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa"
            gtf: "/data/resources/reference/GRCm39/Mus_musculus.GRCm39.113.primary_assembly.gtf"
            version: "mm10"
        ...
    ```
    
    You can download the corresponding files from databases such as NCBI or UCSC, and then edit the relevant entries in the config.yaml file.

- The folder path for the reference files required by JAFFAL when performing gene fusion detection.
    
    ```yaml
    jaffal_refBase: "/data/resources/reference/JAFFA_REFERENCE_FILES_HG38_GENCODE22.V2/"
    ```
    
    You can directly download pre-built reference files provided by the software (currently available only for human hg38/hg19 and mouse mm10) from:https://github.com/Oshlack/JAFFA/wiki/Download#reference-files

    Alternatively, you can follow the official tutorial to build the reference files from scratch: https://github.com/Oshlack/JAFFA/wiki/FAQandTroubleshooting#how-can-i-generate-the-reference-files-for-a-non-supported-genome

    Note: The three analysis modules provided in this pipeline are independent of each other and can be executed selectively. If you do not plan to perform gene fusion detection, you can ignore this file.

### Configuration parameters

Configuration parameters are specified and modified in the `config.yaml` file:

- `glycine_opts`: Parameters for Glycine to identify full-length reads. The default is `-5 AAGCAGTGGTATCAACGCAGAGTACATGGG -3 AAGCAGTGGTATCAACGCAGAGTAC -e 0.25,0.4 -s 100,100 -L 0 -u 10 -l 10`.

- `minimap_genome_index_opts`: Parameters for minimap2 to build the reference genome index, suitable for spliced alignment of transcript sequencing reads. The default is `-x splice -k14`.

- `minimap2_opts`: Parameters for minimap2 alignment. Different parameters are required for genome and transcriptome alignment. The default genome alignment parameters are `-ax splice -uf -k14 --secondary=no`, and the default transcriptome alignment parameters are `-ax map-ont --secondary=no`.

- `modules`: Set the status of the three analysis modules to True (enabled) or False (disabled).
    
    ```
    modules:
        module1: True
        module2: True
        module3: True
    ```

### Prepare for sample sheet

sample sheet is a tab-delimited file (read_manifest.tsv) that stores the information needed to set up this pipeline.

Create a sample sheet file based on the actual sample details. The file must include the following four columns:

- 1st column: Run ID

- 2nd column: Sample name. This column will be used to label samples in the analysis results. Note: Avoid using spaces, hyphens, or other illegal characters in the naming, and do not start the name with a number.

- 3rd column: Sample group. Assign control and experimental groups as control and treated, respectively.

- 4th column: Path to the sequencing data. It is recommended to use an absolute path.

Example:

```
2407C01889011   control_1  control 2407C01889011.fastq.gz
2407C02530011   control_2  control 2407C02530011.fastq.gz
2407C02350011   control_3  control 2407C02350011.fastq.gz
2407C02217011   treated_1  treated 2407C02217011.fastq.gz
2407B02542021   treated_2  treated 2407B02542021.fastq.gz
2407C01579021   treated_3  treated 2407C01579021.fastq.gz
```
    
### Output
    
Data Preprocessing and Quality Control Output：

- fl_reads/<sample_name>/*fq.gz

    Identification results of full-length reads, including quality control statistics for both raw sequencing data and full-length reads, such as length and quality metrics.
        

Output of module 1：

- genome_index/

    Reference genome index built by minimap2
        
- genome_alignments/<sample_name>/*bam

    Alignment output results of each sample's data mapped to the reference genome using minimap2
  
- known_transcripts_depth/

    The Profile.png file contains the read coverage depth analysis results for the coding regions (excluding introns) of known transcripts and their upstream/downstream 2 kb regions.
 
- isoquant/OUT/*transcript_models.gtf
    
    Non-redundant transcript structures obtained from transcriptome reconstruction by IsoQuant
    
- sqanti_qc/
    
    The output results of transcriptome reconstruction annotation and quality control (QC) performed using SQANTI3, where the OUT.transcript_models_classification.txt is a key file. This is a tab-delimited text file, with each row representing a transcript and each column containing the annotation and QC metrics provided by SQANTI3. For detailed file format descriptions, please refer to the [official documentation](https://github.com/ConesaLab/SQANTI3/wiki/Understanding-the-output-of-SQANTI3-QC#classifcols).


Output of module 2：

- transcriptome_index/{transcripts.fa, transcriptome_index.mmi}

    Reference transcriptome index built by minimap2. The reference transcriptome is derived from the provided genome annotation file for the corresponding species.
   
- transcriptome_alignments/<sample_name>/*bam
    
    Alignment output results of each sample's data mapped to the reference transcriptome using minimap2
        
- known_transcripts_coverage/<sample_name>/coverage.tsv
    
    Coverage statistics of each sample's data on known transcripts and genes, with genes categorized into protein_coding and lncRNA for separate analysis. The file format is as follows:
    
    ```
    gene_biotype   coverage          feature
    lncRNA         0.140303894436669 transcript
    protein_coding 0.440216948521282 transcript
    lncRNA         0.331156556109297 gene
    protein_coding 0.768564602826996 gene
    ```
        
- count_reads/<sample_name>/quant.sf
    
    Transcript quantification results for each sample obtained using Salmon.
     
- merge_counts/
    
    A merged expression quantification matrix combining transcript quantification results from all samples, with quantification units in read counts (all_counts.tsv) and TPM (all_TPM.tsv).
    
- diff_exp/
    
    Results from differential expression analysis performed using edgeR, including DGE (Differential Gene Expression) and DTE (Differential Transcript Expression).

Output of module 3：
    
- jaffal_fusion/<sample_name>/
    
    Fusion gene identification results for the corresponding samples. The folder contains two files: jaffa_results.csv and jaffa_results.fasta. For a description of the file formats, please refer to the [official documentation](https://github.com/Oshlack/JAFFA/wiki/OutputDescription).

    The file jaffal_fusion/jaffa_results.csv summarizes the fusion gene identification results for all samples.
        
### Usage
    
The analysis of this pipeline needs to be performed in the conda environment created earlier.

```shell
conda activate pipeline_fulllong_transcriptome
```

The specific usage of the pipeline is as follows:

```shell
snakemake \
--snakefile <path to Snakefile> \ # Specify the path to the Snakefile
--configfile <path to config.yaml> \ # Specify the path to the config.yaml
--config \
manifest=<path to manifest.tsv> \# sample sheet file
specie=<specie> \                # The species of the current data should match the species name specified in config.yaml
with_ERCC=<True | False> \       # Is this an ERCC spike-in experiment?
--directory <output directory> \ # set output directort
--use-conda \
-j <num_cores>                   # set the maximum cpu
```

## Pipeline Overview

### 1. Preprocessing of Raw Sequencing Data (Passed Reads) and Full-Length Reads Identification

The self-developed Glycine tool is used for full-length reads identification and splitting. Glycine identifies full-length reads based on the sequence composition and relative structure of dual-end amplification primers (TSO & TRP) and polyA/T, and provides statistical results for full-length cDNA reads. For chimeric reads, Glycine attempts to split out full-length reads that meet the sequence structure requirements. Only full-length reads are used for downstream analysis.

Additionally, this section performs quality control (QC) statistics on the raw sequencing data (passed reads) and full-length reads data, including length, quality, and GC content.

### 2. Module 1: Transcriptome Reconstruction and Novel Isoform Identification

(1) Mapping Reads to the Reference Genome

The reference genome index is built using minimap2, and reads from each sample are mapped to the reference genome. The alignment results are processed using samtools to generate sorted and indexed BAM files. samtools is also used to generate alignment statistics.

(2) Read Coverage Depth and Uniformity Assessment for Known Transcript Coding Regions

Using the provided reference genome annotation file (GTF file), the coding regions of all known transcripts are obtained. deeptools is then used to assess the read coverage depth and uniformity across these coding regions.

(3) Transcriptome Reconstruction, Annotation, and Quality Control

IsoQuant is used to perform reference genome-guided transcriptome reconstruction. Alignment results from all samples are provided to IsoQuant to generate a comprehensive transcriptome reconstruction result integrating information from all samples.

The transcriptome reconstruction results are then passed to SQANTI3 for transcript structure annotation and quality control. By comparing the reconstructed transcripts with the reference genome annotation, SQANTI3 determines the completeness of each reconstructed transcript relative to known transcript structures, identifies novel isoforms, and evaluates the quality of transcript structures.

### 3. Module 2: Quantification of Known Transcripts and Differential Expression/Splicing Analysis

(1) Construction of the Reference Transcriptome

Using gffread, the mRNA sequences are extracted from the provided reference genome annotation (GTF) and reference genome sequence to construct the reference transcriptome.

(2) Mapping Reads to the Reference Transcriptome

The reference transcriptome index is built using minimap2, and reads from each sample are mapped to the reference transcriptome. The alignment results are processed using samtools to generate sorted and indexed BAM files. samtools is also used to generate alignment statistics.

(3) Statistics on Gene Coverage and Transcript Coverage for Known Genes (Protein-Coding and lncRNA Genes Only)

Using mosdepth, coverage statistics for each transcript are obtained based on the alignment results of each sample to the reference transcriptome. Statistics on gene coverage and transcript coverage are calculated for known genes (protein-coding and lncRNA genes only).

Transcript coverage: The proportion of covered transcripts among all known transcripts.

Gene coverage: The proportion of covered genes among all known genes.

A transcript is considered covered if at least 80% of its sequence is covered. A gene is considered covered if at least one of its transcripts is covered.

(4) Transcript Quantification

Salmon is used to quantify transcripts, generating quantification results for each transcript in each sample. The quantification information from all samples is merged to create a transcript quantification matrix.

(5) Differential Expression Analysis

edgeR is used to perform differential expression analysis, including differential gene expression (DGE) and differential transcript expression (DTE).

(6) Differential Splicing Analysis

SUPPA2 is used to perform differential splicing analysis. The analysis is conducted at two levels:

Local alternative splicing events: These refer to specific splicing patterns within localized regions of a gene, such as skipped exons (SE), retained introns (RI), etc., involving only a small segment of the gene (e.g., a single exon or intron).

Transcript-level events: These focus on expression changes of entire transcripts (isoforms).

### 4. Module 3: Fusion Gene Detection

JAFFA is used to detect fusion genes. The fusion gene detection analysis is performed independently for each sample. After the analysis, the results from all samples are summarized.