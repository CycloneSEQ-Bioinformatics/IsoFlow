#!/usr/bin/env Rscript

suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(tidyr))
suppressMessages(library(pheatmap))
suppressMessages(library(pafr))
suppressMessages(library(GenomicRanges))
suppressMessages(library(Gviz))

args <- commandArgs(trailingOnly=TRUE)

jaffa_results <- args[1]
min_span_reads <- as.numeric(args[2])
ref_version <- args[3]
outdir <- args[4]

if(ref_version %in% c('hg19', 'hg38', 'mm9', 'mm10')){
    if(ref_version == "hg19"){
        suppressMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))
        txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    }else if(ref_version == "hg38"){
        suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
        txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    }else if(ref_version == "mm9"){
        suppressMessages(library(TxDb.Mmusculus.UCSC.mm9.knownGene))
        txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
    }else{
        suppressMessages(library(TxDb.Mmusculus.UCSC.mm10.knownGene))
        txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
    }
}else{
    stop("Only allowed to set ref_version in hg19/hg38/mm9/mm10")
}

inTab <- read.csv(jaffa_results, stringsAsFactors=FALSE)
inTab <- inTab %>% filter(classification=="HighConfidence") %>% rowwise %>% mutate(sample=str_split_1(sample, "\\.")[1])  %>% as.data.frame
inTab_mat <- inTab %>% dplyr::select(fusion.genes, sample, spanning.reads) %>% group_by(fusion.genes, sample) %>% summarise(spanning.reads=sum(spanning.reads)) %>% spread(sample, spanning.reads, fill=0) %>% tibble::column_to_rownames("fusion.genes")
inTab_mat <- inTab_mat[rowMeans(inTab_mat)>1,]
inTab_mat_gather <- inTab_mat %>% tibble::rownames_to_column("fusion") %>% gather(key="sample", value="count", -fusion) %>% mutate(logCount=log(count+0.1))    


print("plot heatmap for spanning reads count")
png(str_glue('{outdir}/span_reads_heatmap.png'))
# pheatmap(inTab_mat, scale='none', cluster_rows=FALSE, cluster_cols=FALSE, display_numbers=TRUE, number_format="%d", main="Spanning reads count (high confidence)")
ggplot(inTab_mat_gather, aes(x=sample, y=fusion, fill=logCount)) + geom_tile() + geom_text(aes(label = round(count, 1)), color = "black", size = 3) + scale_fill_gradient(low = "white", high = "red") + labs(x="", y="", title="Spanning reads count (high confidence)") + theme(legend.position='none')
dev.off()

alignments_gr_list <- list()
# load paf file for each sample
for(sample_name in colnames(inTab_mat)){
    print(str_glue("load paf file for {sample_name}"))
    paf_path <- str_glue("jaffal_fusion/{sample_name}/{sample_name}.usable.fq/{sample_name}.usable.fq_genome.paf")
    paf_data <- read_paf(paf_path) %>% as.data.frame
    alignments_gr <- GRanges(
        seqnames = Rle(paf_data$tname),
        ranges = IRanges(start = paf_data$tstart + 1, width = paf_data$alen),
	strand = Rle(paf_data$strand),
        name = paf_data$qname
    )
    alignments_gr_list[[sample_name]] <- alignments_gr
}

plot_fusionTrack <- function(gene1, chrom1, pos1, gene2, chrom2, pos2, ref_version, alignments_gr, txdb, outfile){
    fusion.genes <- str_glue("{gene1}:{gene2}")
    fusion_breakpoints <- data.frame(chrom1 = chrom1, start1 = pos1-10, end1 = pos1+10, chrom2 = chrom2, start2 = pos2-10, end2 = pos2+10, name = fusion.genes)
    gene_track_chr1 <- GeneRegionTrack(
        txdb,
        genome = ref_version,
        chromosome = chrom1,
        name = gene1
        )

    gene_track_chr2 <- GeneRegionTrack(
        txdb,
        genome = ref_version,
        chromosome = chrom2,
        name = gene2
        )

    fusion_track_chr1 <- AnnotationTrack(
        start = fusion_breakpoints$start1,
        end = fusion_breakpoints$end1,
        chromosome = fusion_breakpoints$chrom1,
        genome = ref_version,
        name = gene1,
        shape = "arrow",
        fill = "red"
    )

    fusion_track_chr2 <- AnnotationTrack(
        start = fusion_breakpoints$start2,
        end = fusion_breakpoints$end2,
        chromosome = fusion_breakpoints$chrom2,
        genome = ref_version,
        name = gene2,
        shape = "arrow",
        fill = "red"
    )

    window_size <- 1e4
    chromosome1 <- chrom1
    start1 <- floor(pos1/window_size)*window_size
    end1 <- start1 + window_size

    chromosome2 <- chrom2
    start2 <- floor(pos2/window_size)*window_size
    end2 <- start2 + window_size

    # 1. 创建 AlignmentsTrack 用于展示测序读长
    alignment_track_chr1 <- AlignmentsTrack(
        alignments_gr,
        isPaired = FALSE,  # 如果是双端测序数据
        name = str_glue("Reads ({chrom1})"),
        chromosome = chrom1,  # 染色体名称
        start = start1,      # 起始位置
        end = end1,          # 结束位置
    )

    alignment_track_chr2 <- AlignmentsTrack(
        alignments_gr,
        isPaired = FALSE,  # 如果是双端测序数据
        name = str_glue("Reads ({chrom2})"),
        chromosome = chrom2,  # 染色体名称
        start = start2,      # 起始位置
        end = end2,          # 结束位置
    )

    genomeAxisTrack_chr1 <- GenomeAxisTrack()
    genomeAxisTrack_chr2 <- GenomeAxisTrack()

    library(grid)

    # 定义图片的宽度和高度
    plot_width <- 6  # 图片宽度（单位：英寸）
    plot_height <- 8  # 图片高度（单位：英寸）
    # 动态计算标题大小
    cex_main <- min(plot_width, plot_height) * 0.15  # 根据图片尺寸调整标题大小

    png(outfile, width=plot_width, height=plot_height, units='in', res=300)
    # 4. 绘制整合的基因结构图
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(2, 1)))  # 2 行 1 列布局

    # 绘制 chr1 的基因结构图
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
    plotTracks(
        list(genomeAxisTrack_chr1, gene_track_chr1, fusion_track_chr1, alignment_track_chr1),
        from = start1,
   	to = end1,
	sizes = c(0.7, 3, 0.5, 5),
        chromosome = chromosome1,
        main = str_glue("Gene Structure with Fusion Points and Read Support ({chrom1})"),
        cex.main = cex_main,  # 使用动态计算的标题大小
        add = TRUE
    )
    popViewport()

    # 绘制 chr2 的基因结构图
    pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
    plotTracks(
        list(genomeAxisTrack_chr2, gene_track_chr2, fusion_track_chr2, alignment_track_chr2),
        from = start2,
        to = end2,
	sizes = c(0.7, 3, 0.5, 5),
        chromosome = chromosome2,
        main = str_glue("Gene Structure with Fusion Points and Read Support ({chrom2})"),
        cex.main = cex_main,  # 使用动态计算的标题大小
        add = TRUE
    )

    popViewport()
    dev.off()
}

passed_inTab <- data.frame()
# plot gene fusion alignments
for(fusion in row.names(inTab_mat)){
    for(sample_name in colnames(inTab_mat)){
        # try to plot for this (fusion, sample_name) tuple
        if(inTab_mat[fusion, sample_name]>=min_span_reads){
            cur_inTab <- inTab %>% filter(fusion.genes==fusion & sample==sample_name)
            passed_inTab <- rbind(passed_inTab, cur_inTab)
            for(i in 1:nrow(cur_inTab)){
                if(cur_inTab$spanning.reads[i]>=min_span_reads){
                    fusion.genes <- cur_inTab$fusion.genes[i]
	            gene1 <- str_split_1(fusion.genes, ':')[1]
                    gene2 <- str_split_1(fusion.genes, ':')[2]
                    chrom1 <- cur_inTab$chrom1[i]
                    pos1 <- cur_inTab$base1[i]
                    chrom2 <- cur_inTab$chrom2[i]
                    pos2 <- cur_inTab$base2[i]
                    print(str_glue("Fusion genes for {sample_name} -- {gene1}({chrom1},{pos1}):{gene2}({chrom2},{pos2}), with {cur_inTab$spanning.reads[i]} spanning reads"))

		    outfile <- str_glue("{outdir}/{sample_name}/{gene1}_{chrom1}_{pos1}.{gene2}_{chrom2}_{pos2}.{cur_inTab$spanning.reads[i]}.png")
                    plot_fusionTrack(gene1, chrom1, pos1, gene2, chrom2, pos2, ref_version, alignments_gr_list[[sample_name]], txdb, outfile)
                }
            }
        }
    }
}

max_row <- passed_inTab[which.max(passed_inTab$spanning.reads),]
max_sample <- max_row$sample
max_fusion.genes <- max_row$fusion.genes
max_gene1 <- str_split_1(max_fusion.genes, ':')[1]
max_gene2 <- str_split_1(max_fusion.genes, ':')[2]
max_chrom1 <- max_row$chrom1
max_pos1 <- max_row$base1
max_chrom2 <- max_row$chrom2
max_pos2 <- max_row$base2
max_spanning.reads <- max_row$spanning.reads
max_path <- str_glue("{max_sample}/{max_gene1}_{max_chrom1}_{max_pos1}.{max_gene2}_{max_chrom2}_{max_pos2}.{max_spanning.reads}.png")
writeLines(max_path, str_glue("{outdir}/which_to_show"))
