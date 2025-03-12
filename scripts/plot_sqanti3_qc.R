#!/usr/bin/env Rscript

suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))

args <- commandArgs(trailingOnly=TRUE)
class_file <- args[1]
outfile <- args[2]

inTab <- read.table(class_file, header=TRUE) %>% select(isoform, structural_category, diff_to_TSS, diff_to_TTS, within_CAGE_peak, polyA_motif_found)

inTab <- inTab %>% rowwise %>% mutate(TSS_support=ifelse((!is.na(diff_to_TSS) & abs(diff_to_TSS)<=50) | within_CAGE_peak, TRUE, FALSE), TTS_support=ifelse((!is.na(diff_to_TTS) & abs(diff_to_TTS)<=50) | polyA_motif_found, TRUE, FALSE))

inTab <- inTab %>% rowwise %>% mutate(type=case_when((TSS_support & TTS_support)~"5'+3'-support", (TSS_support & !TTS_support)~"5'-support only", (!TSS_support & TTS_support)~"3'-support only", TRUE~"No support"))

inTab_summarize <- inTab %>% group_by(structural_category, type) %>% count() %>% as.data.frame %>% filter(structural_category %in% c("full-splice_match", "incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog")) %>% mutate(structural_category_simple=case_match(structural_category, "full-splice_match"~"FSM", "incomplete-splice_match"~"ISM", "novel_in_catalog"~"NIC", "novel_not_in_catalog"~"NNIC", .default = structural_category), type=factor(type, levels=c("5'+3'-support", "5'-support only", "3'-support only", "No support")))

png(outfile, width=4, height=3, units='in', res=300)
ggplot(inTab_summarize, aes(x=structural_category_simple, y=n, fill=type)) + geom_bar(stat="identity") + labs(x="Structural category", y="Count", fill="") + theme_bw() + theme(legend.text = element_text(size = 5), legend.position=c(0.8,0.8), legend.background=element_blank()) + scale_fill_brewer(palette = "Set3")
dev.off()
