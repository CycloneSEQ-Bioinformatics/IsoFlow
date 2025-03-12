#!/usr/bin/env Rscript

suppressMessages(library(dplyr))
suppressMessages(library(stringr))

args <- commandArgs(trailingOnly=TRUE)

trx_infoTab_file <- args[1]
infile <- args[2]
by <- args[3]

txdf <- read.table(trx_infoTab_file, header=TRUE)
inTab <- read.table(infile, header=TRUE)
outfile <- sub("(\\.[^.]+)$", paste0('.add_genename', "\\1"), infile)
if(by=="rowname"){
    raw_colnames <- colnames(inTab)
    inTab <- inTab %>% tibble::rownames_to_column("feature") %>% rowwise %>% mutate(gene_id=str_split_1(feature, ';')[1])
    inTab$gene_name <- txdf$gene_name[match(inTab$gene_id, txdf$gene_id)]

    inTab <- inTab %>% select(all_of(c('feature', 'gene_name', raw_colnames)))
    write.table(inTab, outfile, sep='\t', row.names=FALSE, quote=FALSE)
}else{
    raw_colnames <- colnames(inTab)
    raw_colnames <- raw_colnames[-which(raw_colnames==by)]
    inTab$gene_name <- txdf$gene_name[match(inTab[, by], txdf[, by])]
    inTab <- inTab %>% select(all_of(c(by, 'gene_name', raw_colnames)))
    write.table(inTab, outfile, sep='\t', row.names=FALSE, quote=FALSE)
}
