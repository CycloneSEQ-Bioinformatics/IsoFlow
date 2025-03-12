#!/usr/bin/env Rscript
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))

args <- commandArgs(trailingOnly=TRUE)

mergeTab <- data.frame()
for(i in 1:(length(args)-1)){
	filename_split <- str_split_1(args[i], "/")
	sample_name <- filename_split[length(filename_split)-1]
	inTab <- read.table(args[i], header=TRUE, sep='\t') %>% mutate(sample=sample_name)
	mergeTab <- rbind(mergeTab, inTab)
}

mergeTab$gene_biotype <- factor(mergeTab$gene_biotype, levels=c("protein_coding", "lncRNA", "protein_coding & lncRNA"))
mergeTab <- mergeTab %>% rowwise %>% mutate(feature=str_to_sentence(feature))


write.table(mergeTab, str_glue("{args[length(args)]}.tsv"), row.names=FALSE, quote=FALSE)

png(str_glue("{args[length(args)]}.png"), width=4, height=3, units='in', res=300)
ggplot(mergeTab, aes(x=gene_biotype, y=coverage, fill=sample)) + geom_bar(stat="identity", position="dodge") + facet_wrap(~feature) + labs(x="", y="Coverage", fill="Sample") + theme_bw() + theme(axis.text.x=element_text(size=5, angle=90), legend.text=element_text(size=5)) + scale_fill_brewer(palette="Set3")
dev.off()
