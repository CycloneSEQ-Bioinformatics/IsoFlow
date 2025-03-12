#!/usr/bin/env Rscript

suppressMessages(library("DRIMSeq"))
suppressMessages(library("GenomicFeatures"))
suppressMessages(library("edgeR"))
suppressMessages(library("dplyr"))
suppressMessages(library("ggplot2"))

args <- commandArgs(trailingOnly=TRUE)
de_params_file <-  args[1]
sample_sheet <- args[2]
trx_info_file <- args[3]
cout_tsv <- args[4]
outdir <- args[5]

de_params <- read.csv(de_params_file, sep="\t", stringsAsFactors=FALSE)

cts <- as.matrix(read.csv(cout_tsv, sep="\t", row.names="transcript_id", stringsAsFactors=FALSE))

# Set up sample data frame:
#changed this to sample_id
coldata <- read.csv(sample_sheet, row.names="alias", sep="\t", stringsAsFactors=TRUE)

coldata$sample_id <- rownames(coldata)
# check if control condition exists, sets as reference 
if(!"control" %in% coldata$condition)
  stop("sample_sheet.csv does not contain 'control' 
       condition - unable to set reference.")
coldata$condition <- relevel(coldata$condition, ref = "control")

# a .gff annotation file extension may be gff2(gtf) or gff3 so check in files for use of = in the attribute field
# if '=' present it is gff3 if not it is gtf.
# see https://www.ensembl.org/info/website/upload/gff.html
# and http://gmod.org/wiki/GFF2#Converting_GFF2_to_GFF3
#cat("Checking annotation file type.\n")
#lines <- readLines(file(de_params$ref_annotation[1]), n=10000)
# If transcript_id containing '=' (format eg. transcript_id=xxx)
# annotation type is gff3
#check_file_type <- sum(grepl("transcript_id=", lines))
#if (check_file_type != 0){
#    cat("Annotation file type is gff3.\n")
#    annotation_type <- "gff3"
#} else {
#    # otherwise gtf
#    cat("Annotation file type is gtf.\n")
#    annotation_type <- "gtf"
#}

# Transcript_id versions (eg. ENTXXX.1, eg. ENTXXX.2) represent how many times that transcript reference has been changed 
# during its time in the database.
# Not all annotation files include it as part of the transcript_id - notably Ensembl
# The following handles this.
#cat("Checking annotation file for presence of transcript_id versions.\n")
# Get the first transcript_id from the annotation file by parsing
#lines <- readLines(file(de_params$ref_annotation[1]), n=100000)
# Find transcript_ids in first 1000 lines and check if they contain dot (format eg. ENTXXX.1)
#check_version <- sum(grepl("transcript_id[^;]+\\.", lines))
#if (check_version != 0){
#        # we do not need to strip the count file rows if ref_annotation includes versions
#        cat("Annotation file transcript_ids include versions.\n")
#    } else {
#       # otherwise remove the versions
#        rownames(cts) <- lapply(rownames(cts),  sub, pattern = "\\.\\d+$", replacement = "")
#        cat("Annotation file transcript_ids do not include versions so also strip versions from the counts df.\n")
#    }

#cat("Loading annotation database.\n")
#txdb <- makeTxDbFromGFF(de_params$ref_annotation[1],  format = annotation_type)
#txdf <- select(txdb, keys(txdb,"GENEID"), "TXNAME", "GENEID")
#tab <- table(txdf$GENEID)
#txdf$ntx<- tab[match(txdf$GENEID, names(tab))]
txdf <- read.table(trx_info_file, header=TRUE)

cts <- cts[rownames(cts) %in% txdf$transcript_id, ] # FIXME: filter for transcripts which are in the annotation. Why they are not all there? 

# Reorder transcript/gene database to match input counts:
txdf <- txdf[match(rownames(cts), txdf$transcript_id), ]
rownames(txdf) <- NULL

# Create counts data frame:
counts <- data.frame(gene_id=txdf$gene_id, feature_id=txdf$transcript_id, gene_name=txdf$gene_name, cts)

# output unfiltered version of the counts table now we have paired transcripts with gene ids
write.table(counts, file=stringr::str_glue("{outdir}/unfiltered_transcript_counts_with_genes.tsv"), sep="\t", row.names = FALSE, quote=FALSE)

cat("Filtering counts using DRIMSeq.\n")

d <- dmDSdata(counts=counts, samples=coldata)
trs_cts_unfiltered <- counts(d)

d <- dmFilter(d, min_samps_gene_expr = de_params$min_samps_gene_expr[1], min_samps_feature_expr = de_params$min_samps_feature_expr[1],
        min_gene_expr = de_params$min_gene_expr[1], min_feature_expr = de_params$min_feature_expr[1])
trs_cts <- counts(d)
trs_cts$gene_name <- txdf$gene_name[match(trs_cts$gene_id, txdf$gene_id)]
trs_cts <- cbind(trs_cts %>% select(gene_id, feature_id, gene_name), trs_cts %>% select(-gene_id, -feature_id, -gene_name))
write.table(trs_cts, file=stringr::str_glue("{outdir}/filtered_transcript_counts_with_genes.tsv"), sep="\t", row.names = FALSE, quote=FALSE)

cat("Building model matrix.\n")
design <- model.matrix(~condition, data=DRIMSeq::samples(d))

# Sum transcript counts into gene counts:
cat("Sum transcript counts into gene counts.\n")

gene_cts <- trs_cts_unfiltered %>% dplyr::select(-feature_id)  %>% group_by(gene_id) %>% summarise_all(tibble::lst(sum)) %>% data.frame()
colnames(gene_cts) <- gsub("_sum", "", colnames(gene_cts))
gene_cts$gene_name <- txdf$gene_name[match(gene_cts$gene_id, txdf$gene_id)]
gene_cts <- cbind(gene_cts %>% select(gene_id, gene_name), gene_cts %>% select(-gene_id, -gene_name))
write.table(gene_cts, file=stringr::str_glue("{outdir}/all_gene_counts.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# Output count per million of the gene counts using edgeR CPM
gene_cts <- gene_cts %>% select(-gene_name) %>% tibble::column_to_rownames("gene_id")
cpm_gene_counts <- cpm(gene_cts) %>% as.data.frame %>% tibble::rownames_to_column("gene_id")
cpm_gene_counts$gene_name <- txdf$gene_name[match(cpm_gene_counts$gene_id, txdf$gene_id)]
cpm_gene_counts <- cbind(cpm_gene_counts %>% select(gene_id, gene_name), cpm_gene_counts %>% select(-gene_id, -gene_name))
write.table(cpm_gene_counts, file=stringr::str_glue("{outdir}/cpm_gene_counts.tsv"), sep="\t", quote=FALSE, row.names = FALSE)

# Output count per million of the transcript counts using edgeR CPM
trs_cts_unfiltered <- trs_cts_unfiltered  %>% select(-gene_id) %>% tibble::column_to_rownames("feature_id")
cpm_trs_counts <- cpm(trs_cts_unfiltered) %>% as.data.frame %>% tibble::rownames_to_column("transcript_id")
cpm_trs_counts$gene_name <- txdf$gene_name[match(cpm_trs_counts$transcript_id, txdf$transcript_id)]
cpm_trs_counts <- cbind(cpm_trs_counts %>% select(transcript_id, gene_name), cpm_trs_counts %>% select(-transcript_id, -gene_name))
write.table(cpm_trs_counts, file=stringr::str_glue("{outdir}/cpm_trs_counts.tsv"), sep="\t", quote=FALSE, row.names = FALSE)


# Differential gene expression using edgeR:
cat("Running differential gene expression analysis using edgeR.\n")

y <- DGEList(gene_cts)
y <- calcNormFactors(y)
y <- estimateDisp(y,design)
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit)
edger_res <- topTags(qlf, n=nrow(y), sort.by="PValue")[[1]] %>% tibble::rownames_to_column("gene_id")
edger_res$gene_name <- txdf$gene_name[match(edger_res$gene_id, txdf$gene_id)]
edger_res <- cbind(edger_res %>% select(gene_id, gene_name), edger_res %>% select(-gene_id, -gene_name))
write.table(edger_res, file=stringr::str_glue("{outdir}/results_dge.tsv"), sep="\t", row.names = FALSE, quote=FALSE)

# plot MA figure
png(stringr::str_glue("{outdir}/results_dge.maplot.png"), width=4, height=4, units='in', res=300)
# create status vector
qlf$table$label <- ifelse(
  qlf$table$PValue<0.05 & qlf$table$logFC>=1, 
  'Up', 
  ifelse(
    qlf$table$PValue<0.05 & qlf$table$logFC<=-1,
    'Down',
    'No'
  )
)
ggplot(qlf$table, aes(x=logCPM, y=logFC, color=label)) + geom_point(size=0.5) + geom_hline(yintercept=c(-1,1), color="blue", linetype = "dashed") + scale_color_manual(values=c("blue", "gray", "red")) + labs(x="Average log CPM", y="Log fold-change", title="Control vs. Treatment (DGE)", color="") + theme_bw() + theme(legend.position=c(0.8,0.8), legend.background=element_blank(), plot.title=element_text(hjust=0.5))
#plotMD(qlf, status=qlf$table$label, values=c("Up","Down","No"), hl.col=c("red","blue","black"), main="Control vs. Treatment (DGE)")
#abline(h=c(-1,1), col="blue", lty = "dashed")
#plotQLDisp(fit)
dev.off()

# plot volcano figure
png(stringr::str_glue("{outdir}/results_dge.volcanoplot.png"), width=4, height=4, units='in', res=300)
ggplot(qlf$table, aes(x=logFC, y=-log10(PValue), color=label)) + geom_point(size=0.5) + scale_color_manual(values=c("blue", "gray", "red")) + labs(title="Control vs. Treatment (DGE)", color="") + theme_bw() + theme(legend.position=c(0.8,0.8), legend.background=element_blank(), plot.title=element_text(hjust=0.5))
dev.off()

cat("Running differential transcript expression analysis using edgeR.\n")
trs_cts <- trs_cts %>% tibble::column_to_rownames("feature_id") %>% select(-gene_id)
y <- DGEList(trs_cts)
y <- calcNormFactors(y)
y <- estimateDisp(y,design)
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit)
edger_res <- topTags(qlf, n=nrow(y), sort.by="PValue")[[1]] %>% tibble::rownames_to_column("transcript_id")
edger_res$gene_name <- txdf$gene_name[match(edger_res$transcript_id, txdf$transcript_id)]
edger_res <- cbind(edger_res %>% select(transcript_id, gene_name), edger_res %>% select(-transcript_id, -gene_name))
write.table(edger_res, file=stringr::str_glue("{outdir}/results_dte.tsv"), sep="\t", row.names = FALSE, quote=FALSE)

# plot MA figure
png(stringr::str_glue("{outdir}/results_dte.maplot.png"), width=4, height=4, units='in', res=300)
# create status vector
qlf$table$label <- ifelse(
  qlf$table$PValue<0.05 & qlf$table$logFC>=1,
  'Up',
  ifelse(
    qlf$table$PValue<0.05 & qlf$table$logFC<=-1,
    'Down',
    'No'
  )
)
ggplot(qlf$table, aes(x=logCPM, y=logFC, color=label)) + geom_point(size=0.5) + geom_hline(yintercept=c(-1,1), color="blue", linetype = "dashed") + scale_color_manual(values=c("blue", "gray", "red")) + labs(x="Average log CPM", y="Log fold-change", title="Control vs. Treatment (DTE)", color="") + theme_bw() + theme(legend.position=c(0.8,0.8), legend.background=element_blank(), plot.title=element_text(hjust=0.5))
#plotMD(qlf, status=qlf$table$label, values=c("Up","Down","No"), hl.col=c("red","blue","black"), main="Control vs. Treatment (DTE)")
#abline(h=c(-1,1), col="blue", lty = "dashed")
dev.off()

# plot volcano figure
png(stringr::str_glue("{outdir}/results_dte.volcanoplot.png"), width=4, height=4, units='in', res=300)
ggplot(qlf$table, aes(x=logFC, y=-log10(PValue), color=label)) + geom_point(size=0.5) + scale_color_manual(values=c("blue", "gray", "red")) + labs(title="Control vs. Treatment (DTE)", color="") + theme_bw() + theme(legend.position=c(0.8,0.8), legend.background=element_blank(), plot.title=element_text(hjust=0.5))
dev.off()

# Differential transcript usage using DEXSeq:
#suppressMessages(library("DEXSeq"))
#cat("Running differential transcript usage analysis using DEXSeq.\n")

#sample.data<-DRIMSeq::samples(d)
#count.data <- round(as.matrix(counts(d)[,-c(1:2)]))
#dxd <- DEXSeqDataSet(countData=count.data, sampleData=sample.data, design=~sample + exon + condition:exon, featureID=trs_cts$feature_id, groupID=trs_cts$gene_id)
#dxd <- estimateSizeFactors(dxd)
#dxd <- estimateDispersions(dxd)
#dxd <- testForDEU(dxd, reducedModel=~sample + exon)
#dxd <- estimateExonFoldChanges( dxd, fitExpToVar="condition")
#dxr <- DEXSeqResults(dxd, independentFiltering=FALSE)

#dev.off()
#pdf(stringr::str_glue("{outdir}/results_dtu.pdf"))
#plotMA(dxr, cex=0.8, alpha=0.05) 
#plotDispEsts(dxd)
#dev.off()
#
#qval <- perGeneQValue(dxr) 
#dxr.g<-data.frame(gene=names(qval), qval)
#dxr.g <- dxr.g[order(dxr.g$qval),]
#
#dxr_out <- as.data.frame(dxr[,c("featureID", "groupID", "pvalue")])
#dxr_out <- dxr_out[order(dxr$pvalue),]
#
#write.table(dxr.g, file=stringr::str_glue("{outdir}/results_dtu_dexseq_gene.tsv"), row.names=FALSE, sep="\t", quote=FALSE)
#write.table(dxr_out, file=stringr::str_glue("{outdir}/results_dtu_dexseq_transcript.tsv"), row.names=FALSE, sep="\t", quote=FALSE)
#
## and writing out some of the DEXSeq metrics
#colnames(dxr)[grep("log2fold", colnames(dxr))] <- "log2fold"
#MADTUdata <- data.frame(dxr)[order(dxr$padj),c("exonBaseMean", "log2fold", "pvalue", "padj")]
#MADTUdata$exonBaseMean <- log2(MADTUdata$exonBaseMean)
#colnames(MADTUdata)[which(colnames(MADTUdata)=="exonBaseMean")] <- "Log2MeanExon"
#colnames(MADTUdata)[which(colnames(MADTUdata)=="log2fold")] <- "Log2FC"
#write.table(MADTUdata %>% tibble::rownames_to_column("featureID"), file=stringr::str_glue("{outdir}/results_dtu_dexseq_transcript.detail.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
#
## stageR analysis of DEXSeq results:
#cat("stageR analysis\n")
#library(stageR)
#
#cat("Running stageR analysis on the differential transcript usage results.\n")
#pConfirmation <- matrix(dxr$pvalue, ncol=1)
#
#dimnames(pConfirmation) <- list(dxr$featureID, "transcript")
#pScreen <- qval
#tx2gene <- as.data.frame(dxr[,c("featureID", "groupID")])
#
#stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation, pScreenAdjusted=TRUE, tx2gene=tx2gene)
## note: the choice of 0.05 here means you can *only* threshold at 5% OFDR later
#stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.10)
#suppressWarnings({dex.padj <- getAdjustedPValues(stageRObj, order=FALSE, onlySignificantGenes=FALSE)})
#
## dex.padj <- dex.padj[,-1]
#write.table(dex.padj %>% rename(p_gene=gene, p_transcript=transcript), file=stringr::str_glue("{outdir}/results_dtu_stageR.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
