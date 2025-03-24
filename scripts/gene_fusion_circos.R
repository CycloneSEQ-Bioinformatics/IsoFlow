#!/usr/bin/env Rscript
suppressMessages(library(RCircos))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))

args <- commandArgs(trailingOnly = TRUE)
cyto_band_ideo_path <- args[1]
gene_info_path <- args[2]
fusion_info_path <- args[3]
outdir <- args[4]

cyto_band_ideo <- read.csv(cyto_band_ideo_path, sep = "\t")
gene_info <- read.csv(gene_info_path, sep = "\t")
fusion_info <- read.csv(fusion_info_path, sep = ",")

fusion_info <- fusion_info %>% 
    filter(classification == "HighConfidence") %>% 
    separate(col = fusion.genes, into = c("gene1", "gene2"), sep = ":")

fusion_info <- fusion_info %>% 
    left_join(gene_info, by = c("gene1" = "Gene"), suffix = c("", "_gene1")) %>% 
    left_join(gene_info, by = c("gene2" = "Gene"), suffix = c("_gene1", "_gene2")) %>%
    mutate(
        chromStart_gene1 = ifelse(is.na(chromStart_gene1), base1, chromStart_gene1),
        chromEnd_gene1 = ifelse(is.na(chromEnd_gene1), base1, chromEnd_gene1),
        chromStart_gene2 = ifelse(is.na(chromStart_gene2), base2, chromStart_gene2),
        chromEnd_gene2 = ifelse(is.na(chromEnd_gene2), base2, chromEnd_gene2),
    )

fusion_info <- fusion_info %>%
    mutate(
        chromStart_gene1_fusion = ifelse(strand1 == "+", base1, chromStart_gene1),
        chromEnd_gene1_fusion = ifelse(strand1 == "+", chromEnd_gene1, base1),
        chromStart_gene2_fusion = ifelse(strand2 == "+", base2, chromStart_gene2),
        chromEnd_gene2_fusion = ifelse(strand2 == "+", chromEnd_gene2, base2),
    ) %>%
    filter(
        (chromEnd_gene1_fusion - chromStart_gene1_fusion >= 0) & 
        (chromEnd_gene2_fusion - chromStart_gene2_fusion >= 0)
    )

fusion_info <- fusion_info %>%
    mutate(
        line_width = case_when(
            fusion_info$spanning.reads <= 100 ~ 1,
            fusion_info$spanning.reads > 100 & fusion_info$spanning.reads <= 200 ~ 2,
            fusion_info$spanning.reads > 200 & fusion_info$spanning.reads <= 300 ~ 3,
            fusion_info$spanning.reads > 300 & fusion_info$spanning.reads <= 400 ~ 4,
            fusion_info$spanning.reads > 400 & fusion_info$spanning.reads <= 500 ~ 5,
            fusion_info$spanning.reads > 500 ~ 6
        )
    )

fusion_gene1_info <- fusion_info %>%
    select(
        Chromosome = chrom1,
        chromStart = chromStart_gene1,
        chromEnd = chromEnd_gene1,
        Gene = gene1,
    )

fusion_gene2_info <- fusion_info %>%
    select(
        Chromosome = chrom2,
        chromStart = chromStart_gene2,
        chromEnd = chromEnd_gene2,
        Gene = gene2,
    )

fusion_gene_label <- bind_rows(fusion_gene1_info, fusion_gene2_info) %>% distinct()

fusion_gene_plot <- fusion_info %>%
    select(
        Chromosome = chrom1,
        chromStart = chromStart_gene1_fusion,
        chromEnd = chromEnd_gene1_fusion,
        Chromosome.1 = chrom2,
        chromStart.1 = chromStart_gene2_fusion,
        chromEnd.1 = chromEnd_gene2_fusion,
    ) %>%
    mutate(
        PlotColor = case_when(
            fusion_info$known == "Yes" ~ "blue",
	    fusion_info$known == "-" ~ "red",
        )
    )

fusion_gene_line_width <- fusion_info[["line_width"]]

png(str_glue("{outdir}/gene_fusion_circos_plot.png"), width = 4200, height = 4200, res = 600)

RCircos.Set.Core.Components(cyto_band_ideo, chr.exclude = NULL, tracks.inside = 20, tracks.outside = 0)
RCircos.Set.Plot.Area(margins = 0)
RCircos.Chromosome.Ideogram.Plot()

RCircos.Gene.Connector.Plot(fusion_gene_label, track.num = 1, side = "in")
RCircos.Gene.Name.Plot(fusion_gene_label, name.col = 4, track.num = 2, side = "in")

RCircos.Link.Plot(fusion_gene_plot, track.num = 3, lineWidth = fusion_gene_line_width)

dev.off()
