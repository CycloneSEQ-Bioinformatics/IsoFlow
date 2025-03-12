# Bulk全长转录组分析流程使用说明文档

基本介绍

![](pictures/IsoFlow_overview.png)

该流程用于CycloneSEQ bulk 全长转录组分析，分析流程分为三个独立的模块：

**module1**：转录组重构与新型转录本鉴定（transcriptome reconstruction and novel isoform discovery）

**module2**：已知转录本的定量与差异表达/剪切分析（transcript quantification and differential expression/splicing analysis）

**module3**：融合基因检测（gene fusion detection）

默认执行这三个模块的分析任务，若用户只想执行其中的部分模块，可以通过修改配置文件`config.yaml`，进行对应模块的开启或关闭

安装与使用

- 安装
    
    （1）安装miniconda，参考 [官方指导](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) 完成安装
    
    （2）将本流程的文件夹，拷贝到本地
    
    （3）创建用于该分析流程的conda环境
    
    ```shell
    conda env create -n pipeline_fulllong_transcriptome -f env.yaml
    ```

    （4）SUPPA2的安装

    ```shell
    wget -O SUPPA-2.4.tar.gz https://github.com/comprna/SUPPA/archive/refs/tags/v2.4.tar.gz
    tar zxvf SUPPA-2.4.tar.gz
    ````

    并修改config.yaml，以提供SUPPA的安装路径

    ```shell
    SUPPA_HOME: "/PATH/TO/SUPPA-2.4"
    ```
    
    （5）JAFFA的安装
    
    JAFFA为用于融合基因检测的工具，其运行环境相对独特，需要为其专门创建conda环境，并在该环境中完成工具的安装，具体流程如下：
    
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
    
    （6）SQANT3的安装
    
    SQANT3为用于转录组注释与指控的工具，其运行环境相对独特，需要为其专门创建conda环境，并在该环境中完成工具的安装，具体流程如下：
    
    ```bash
    wget -O SQANTI3-5.2.2.tar.gz https://github.com/ConesaLab/SQANTI3/archive/refs/tags/v5.2.2.tar.gz
    tar zxvf SQANTI3-5.2.2.tar.gz
    cd SQANTI3-5.2.2
    
    conda env create -f SQANTI3.conda_env.yml
    
    # 查看和测试该conda环境是否创建成功
    conda env list # 查看是否存在SQANTI3.env的环境
    conda activate SQANTI3.env # 激活SQANTI3.env的环境
    ```

    （7）ggsashimi的安装

    ggsashimi为用于差异可变剪切可视化的工具，其运行环境相对独特，需要为其专门创建conda环境，并在该环境中完成工具的安装，具体流程如下：

    ```shell
    conda config --set channel_priority flexible
    conda env create -n ggsashimi_env -f env.ggsashimi.yaml
    ```
    
    （8）Glycine的安装
    
    Glycine为CycloneSEQ自研的全长reads鉴定工具，可以从下机测序数据中，依据双端扩增引物（TSO & TRP）和polyA/T的序列组成及其相对结构，识别全长reads
    
    由于获得的软件为可执行文件，只需将其保存到某个文件下，并将其加入环境变量即可

    ```shell
    wget -c https://github.com/CycloneSEQ-Bioinformatics/Glycine/releases/download/v1.0.0/glycine.tar.gz
    tar zxvf glycine.tar.gz
    ```

    
    
- 输入
    
    通用的输入文件，例如参考基因组等，和配置参数，由`config.yaml`文件指定：
    
    输入文件的准备与指定
    
    - （多个）目标物种的参考基因组序列和注释
        
        ```jsx
        reference:
            Human:
                fasta: "/data/resources/reference/hg38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
                gtf: "/data/resources/reference/hg38/Homo_sapiens.GRCh38.113.primary_assembly.gtf"
            Mouse:
                fasta: "/data/resources/reference/GRCm39/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa"
                gtf: "/data/resources/reference/GRCm39/Mus_musculus.GRCm39.113.primary_assembly.gtf"
            ...
        ```
        
        可以从NCBI、UCSC等数据库下载对应文件，并编辑`config.yaml`文件的对应项
        
    - JAFFAL进行融合基因检测时，需要用到的参考文件的文件夹路径
        
        ```jsx
        jaffal_refBase: "/data/resources/reference/JAFFA_REFERENCE_FILES_HG38_GENCODE22.V2/"
        ```
        
        可以直接现在软件官方提供的已构建好的参考文件（仅提供人类hg38/hg19，和小鼠mm10下的）：[https://github.com/Oshlack/JAFFA/wiki/Download#reference-files](https://github.com/Oshlack/JAFFA/wiki/Download#reference-files)
        
        也可以按照软件官方提供的教程，自行从头构建参考文件：[https://github.com/Oshlack/JAFFA/wiki/FAQandTroubleshooting#how-can-i-generate-the-reference-files-for-a-non-supported-genome](https://github.com/Oshlack/JAFFA/wiki/FAQandTroubleshooting#how-can-i-generate-the-reference-files-for-a-non-supported-genome)
        
        注意：本流程提供的三个分析模块是相互独立的，可以选择行执行，若不准备进行融合基因检测，可以忽略该文件
        
    
    配置参数设定
    
    - `glycine_opts`：设定Glycine进行全长转录本鉴定的参数，默认为`-5 AAGCAGTGGTATCAACGCAGAGTACATGGG -3 AAGCAGTGGTATCAACGCAGAGTAC -e 0.25,0.4 -s 100,100 -L 0 -u 10 -l 10`
    - `minimap_genome_index_opts`：设定minimap2构建参考基因组索引，以适用于转录本测序片段的剪切跳跃式比对，默认为`-x splice -k14`
    - `minimap2_opts`：设定minimap2进行比对的参数，需要对基因组比对和转录组比对设定不同的参数，默认基因组比对参数为`-ax splice -uf -k14 --secondary=no`，转录组比对参数为`-ax map-ont --secondary=no`
    - modules：设定三个分析模块的开启（True）或关闭（False）状态
        
        ```jsx
        modules:
            module1: True # 转录组重构与新型转录本鉴定
            module2: True # 已知转录本的定量与差异表达分析
            module3: True # 融合基因检测
        ```
        
    
    样本信息文件（`read_manifest.tsv`）
    
    根据实际样本信息，创建样本信息文件，要求包含以下4列：
    
    - 第一列：run id
    - 第二列：样本名，分析结果中将用此列进行样本标记
    - 第三列：样本分组，将对照组与实验组分别设为`control`和`treated`
    - 第四列：测序数据保存路径，建议使用绝对路径
    
    ```
    2407C01889011   control_1  control 2407C01889011.fastq.gz
    2407C02530011   control_2  control 2407C02530011.fastq.gz
    2407C02350011   control_3  control 2407C02350011.fastq.gz
    2407C02217011   treated_1  treated 2407C02217011.fastq.gz
    2407B02542021   treated_2  treated 2407B02542021.fastq.gz
    2407C01579021   treated_3  treated 2407C01579021.fastq.gz
    ```
    
- 输出
    
    数据预处理及指控输出：
    
    - fl_reads/<sample_name>/*fq.gz
        
        全长reads鉴定结果，并包含下机数据与全长reads的长度和质量指控统计结果
        
    
    模块1的输出：
    
    - genome_index/genome_splice_index.mmi
        
        minimap2构建的参考基因组索引
        
    - genome_alignments/<sample_name>/*bam
        
        各样本数据用minimap2比对到参考基因组的比对输出结果
        
    - known_transcripts_depth/
        
        其中的Profile.png，为已知转录本编码区间（不包括intron）及其上下游2kb区间的read覆盖深度分析结果
        
    - isoquant/OUT/*transcript_models.gtf
        
        IsoQuant执行转录组重构得到的非冗余转录本结构
        
    - sqanti_qc/
        
        使用SQANTI3，进行转录组重构结果的注释和QC输出结果
        
    
    模块2的输出：
    
    - transcriptome_index/{transcripts.fa, transcriptome_index.mmi}
        
        minimap2构建的参考转录组索引。参考转录组来源于提供的对应物种的基因组注释
        
    - transcriptome_alignments/<sample_name>/*bam
        
        各样本数据用minimap2比对到参考转录组的比对输出结果
        
    - known_transcripts_coverage/<sample_name>/coverage.tsv
        
        各样本数据对已知转录本和基因的覆盖情况，将基因分成protein_coding和lncRNA分别进行统计，文件格式如下
        
        ```
        gene_biotype   coverage          feature
        lncRNA         0.140303894436669 transcript
        protein_coding 0.440216948521282 transcript
        lncRNA         0.331156556109297 gene
        protein_coding 0.768564602826996 gene
        ```
        
    - count_reads/<sample_name>/quant.sf
        
        各样本的转录本定量结果
        
    - merge_counts/all_counts.tsv
        
        合并所有样本的转录本定量结果，得到的表达定量矩阵，定量单位为read count
        
    - de_analysis/
        
        差异表达分析结果，包含DGE（Differential Gene Expression）和DTU（Differential Transcript Usage）
        
    
    模块3的输出：
    
    - jaffal_fusion/<sample_name>/
        
        融合基因鉴定结果，文件夹下包含两个文件`jaffa_results.csv`和`jaffa_results.fasta`，文件格式说明见 [官方文档](https://github.com/Oshlack/JAFFA/wiki/OutputDescription)
        
- 使用
    
    需要在前面创建的conda环境下进行本流程的分析
    
    ```shell
    conda activate pipeline_fulllong_transcriptome
    ```
    
    流程的具体用法如下：
    
    ```shell
    snakemake \
    --snakefile <path to Snakefile> \ # 指定Snake文件的路径
    --configfile <path to config.yaml> \ # 指定config.yaml文件的路径
    --config \
    manifest=<path to manifest.tsv> \# 样本信息文件
    specie=<specie> \                # 当前数据的物种，请与config.yaml中的物种名一致
    with_ERCC=<True | False> \       # 是否为ERCC spike-in实验
    --directory <output directory> \ # 指定输出目录
    --use-conda \
    -j <num_cores>                   # 指定最大运行CPU数
    ```