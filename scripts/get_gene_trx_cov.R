#!/usr/bin/env Rscript
suppressMessages(library("GenomicFeatures"))
suppressMessages(library("biomaRt"))
suppressMessages(library(dplyr))
suppressMessages(library(rtracklayer))

args <- commandArgs(trailingOnly=TRUE)
trx_infoTab_file <- args[1]
cov_file <- args[2]
outfile <- args[3]

# cat("Checking annotation file type.\n")
#lines <- readLines(file(ref_annotation), n=10000)
# If transcript_id containing '=' (format eg. transcript_id=xxx)
# annotation type is gff3
# check_file_type <- sum(grepl("transcript_id=", lines))
#if (check_file_type != 0){
#    cat("Annotation file type is gff3.\n")
#    annotation_type <- "gff3"
#} else {
#    # otherwise gtf
#    cat("Annotation file type is gtf.\n")
#    annotation_type <- "gtf"
#}

#cat("Loading annotation database.\n")
#txdb <- makeTxDbFromGFF(ref_annotation, format = annotation_type)
#txdf <- biomaRt::select(txdb, keys(txdb,"GENEID"), "TXNAME", "GENEID")
#gtf <- import(ref_annotation)
#gene_info <- mcols(gtf)[, c("gene_id", "gene_biotype")] %>% as.data.frame %>% rename(GENEID="gene_id") %>% distinct()

txdf <- read.table(trx_infoTab_file, header=TRUE) 

#cov_tab <- read.table(cov_file, col.names=c('transcript_id', 'depth', 'coverage')) %>% filter(depth==1)
cov_tab <- read.table(cov_file, col.names=c('TrxID', 'depth', 'coverage')) %>% filter(depth==1)

#txdf <- txdf %>% merge(cov_tab, by="transcript_id", all=TRUE) %>% mutate(depth=1, coverage=ifelse(is.na(coverage), 0, coverage))
#txdf <- txdf %>% merge(gene_info, by="GENEID", all=TRUE)
txdf <- txdf %>% merge(cov_tab, by="TrxID", all=TRUE) %>% mutate(depth=1, coverage=ifelse(is.na(coverage), 0, coverage))

cov_stat <- rbind(
		txdf %>% dplyr::filter(gene_biotype %in% c("protein_coding", "lncRNA")) %>% group_by(gene_biotype) %>% summarise(coverage=mean(coverage>=0.8)) %>% mutate(feature="transcript"),
		txdf %>% dplyr::filter(gene_biotype %in% c("protein_coding", "lncRNA")) %>% group_by(GeneID, gene_biotype) %>% summarise(coverage=max(coverage)) %>% as.data.frame %>% group_by(gene_biotype) %>% summarise(coverage=mean(coverage>=0.8)) %>% mutate(feature="gene"),
		txdf %>% dplyr::filter(gene_biotype %in% c("protein_coding", "lncRNA")) %>% summarise(coverage=mean(coverage>=0.8)) %>% mutate(gene_biotype="protein_coding & lncRNA", feature="transcript"),
		txdf %>% dplyr::filter(gene_biotype %in% c("protein_coding", "lncRNA")) %>% group_by(GeneID) %>% summarise(coverage=max(coverage)) %>% as.data.frame %>% summarise(coverage=mean(coverage>=0.8))  %>% mutate(gene_biotype="protein_coding & lncRNA", feature="gene")
		  )
write.table(cov_stat, outfile, sep='\t', row.names=FALSE, quote=FALSE)
