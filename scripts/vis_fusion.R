library(Gviz)
library(GenomicAlignments)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)

# 获取命令行参数
args <- commandArgs(trailingOnly = TRUE)
paf_file <- args[1]
gtf_file <- args[2]
chrom1 <- args[3]
pos1 <- as.integer(args[4])
gene1 <- args[5]
chrom2 <- args[6]
pos2 <- as.integer(args[7])
gene2 <- args[8]
outfile <- args[9]

paf_data <- read.table(paf_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE, fill = TRUE)

# 2. 提取必要列（假设 PAF 文件的前 12 列是标准格式）
colnames(paf_data)[1:12] <- c(
    "query_name", "query_length", "query_start", "query_end", "strand",
    "target_name", "target_length", "target_start", "target_end",
    "num_matches", "alignment_length", "mapping_quality"
)

# 3. 将 PAF 数据转换为 GRanges 对象
alignments_gr <- GRanges(
    seqnames = Rle(paf_data$target_name),  # 染色体名称
    ranges = IRanges(start = paf_data$target_start + 1, width = paf_data$alignment_length),  # 起始位置和长度
    strand = Rle(paf_data$strand),         # 链方向
    name = paf_data$query_name             # 读长名称
)

gene_track_chr1 <- GeneRegionTrack(
    gtf_file,
    genome = "hg38",  # 参考基因组版本
    chromosome = chrom1,  # 染色体名称
    name = gene1
)

gene_track_chr2 <- GeneRegionTrack(
    gtf_file,
    genome = "hg38",  # 参考基因组版本
    chromosome = chrom2,  # 染色体名称
    name = gene2
)

# 2. 创建 AnnotationTrack 用于标注融合断点
fusion_track_chr1 <- AnnotationTrack(
    start = fusion_breakpoints$start1,
    end = fusion_breakpoints$end1,
    chromosome = fusion_breakpoints$chrom1,
    genome = "hg38",
    name = gene1,
    shape = "arrow",
    fill = "red"
)

fusion_track_chr2 <- AnnotationTrack(
    start = fusion_breakpoints$start2,
    end = fusion_breakpoints$end2,
    chromosome = fusion_breakpoints$chrom2,
    genome = "hg38",
    name = gene2,
    shape = "arrow",
    fill = "red"
)


# 1. 定义基因组区域
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
    stacking = "full",  # 禁用堆叠
    size = 0.5,         # 调整轨道高度
    min.height = 0.1,  # 最小高度
    max.height = 0.8   # 最大高度
)

alignment_track_chr2 <- AlignmentsTrack(
    alignments_gr,
    isPaired = TRUE,  # 如果是双端测序数据
    name = str_glue("Reads ({chrom2})"),
    chromosome = chrom2,  # 染色体名称
    start = start2,      # 起始位置
    end = end2,          # 结束位置
    stacking = "full",  # 禁用堆叠
    size = 0.5,         # 调整轨道高度
    min.height = 0.1,  # 最小高度
    max.height = 0.8   # 最大高度
)

# 2. 创建基因组坐标轴轨道
genomeAxisTrack_chr1 <- GenomeAxisTrack()
genomeAxisTrack_chr2 <- GenomeAxisTrack()

library(grid)

# 定义图片的宽度和高度
plot_width <- 6  # 图片宽度（单位：英寸）
plot_height <- 8  # 图片高度（单位：英寸）

# 动态计算标题大小
cex_main <- min(plot_width, plot_height) * 0.2  # 根据图片尺寸调整标题大小

# 4. 绘制整合的基因结构图
pdf(outfile, width = plot_width, height = plot_height)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1)))  # 2 行 1 列布局

# 绘制 chr1 的基因结构图
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
plotTracks(
    list(genomeAxisTrack_chr1, gene_track_chr1, fusion_track_chr1, alignment_track_chr1),
    from = start1,
    to = end1,
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
    chromosome = chromosome2,
    main = str_glue("Gene Structure with Fusion Points and Read Support ({chrom2})"),
    cex.main = cex_main,  # 使用动态计算的标题大小
    add = TRUE
)
popViewport()
dev.off()
