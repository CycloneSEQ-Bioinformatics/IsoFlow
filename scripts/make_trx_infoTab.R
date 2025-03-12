#!/usr/bin/env Rscript
suppressMessages(library(dplyr))
suppressMessages(library(rtracklayer))

args <- commandArgs(trailingOnly=TRUE)
ref_annotation <- args[1]
outfile <- args[2]

#cat("Checking annotation file type.\n")
#lines <- readLines(file(ref_annotation), n=10000)
## If transcript_id containing '=' (format eg. transcript_id=xxx)
## annotation type is gff3
#check_file_type <- sum(grepl("transcript_id=", lines))
#if (check_file_type != 0){
#    cat("Annotation file type is gff3.\n")
#    annotation_type <- "gff3"
#} else {
#    # otherwise gtf
#    cat("Annotation file type is gtf.\n")
#    annotation_type <- "gtf"
#}
#
#cat("Loading annotation database.\n")
#txdb <- makeTxDbFromGFF(ref_annotation, format = annotation_type)
#txdf <- biomaRt::select(txdb, keys(txdb,"GENEID"), "TXNAME", "GENEID")
gtf <- import(ref_annotation)
gene_info <- mcols(gtf)[, c("type", "gene_id", "transcript_id", "gene_biotype", "gene_name")] %>% as.data.frame %>% filter(type=="transcript") %>% as.data.frame %>% select(-type) %>% distinct()

write.table(gene_info, outfile, sep='\t', row.names=FALSE, quote=FALSE)
