
#
# R script for processing data from proteome mass spec analysis
# after Signal Peptide prediction
#
# orignally by Ryohei Thomas Nakano, PhD
# nakano@mpipz.mpg.de
#

options(warn=-1)

# cleaning up
rm(list=ls())

# libraries
library(stringr,   quietly=T, warn.conflicts=F)
library(reshape2,  quietly=T, warn.conflicts=F)
library(dplyr,     quietly=T, warn.conflicts=F)
library(pipeR,     quietly=T, warn.conflicts=F)


# directories
original_data  <- "/your/data/directory/original_data/"
processed_data <- "/your/data/directory/processed_data/"
fig            <- "/your/data/directory/fig/"
stat           <- "/your/data/directory/statistics/"


# data import
signalP  <- read.table(paste(processed_data, "signalP.txt", sep=""),                  header=T, sep="", comment.char="",stringsAsFactors=F)
LFQ_melt <- read.table(paste(processed_data, "LFQ_melt_table_corrected.txt", sep=""), header=T, stringsAsFactors=F)
design   <- read.table(paste(processed_data, "design.txt", sep=""),                   header=T, stringsAsFactors=F)

# convert LFQ to dcast matrix
idx <- names(LFQ_melt) %in% c("SampleID", "LFQ")
dcast <- dcast(LFQ_melt[,!idx], ID ~ genotype+timepoint+rep, value.var="mLFQ")
LFQ <- as.matrix(dcast[,-1])
colnames(LFQ) <- names(dcast)[-1]
rownames(LFQ) <- dcast$ID


# sort
idx <- order(signalP$name)
signalP <- signalP[idx, ]

# aggregate by genes
signalP$gene <- str_replace(signalP$name, "\\..*", "")
summary <- signalP %>>% group_by(gene) %>>% summarise(count=sum(X. == "Y"))

# genes with at least one variant with predicted SP
idx <- summary$count > 0
SP_genes <- summary$gene[idx]

idx <- signalP$gene %in% SP_genes
signalP_SP <- signalP[idx, ]

idx <- rownames(LFQ) %in% SP_genes

# output
write.table(SP_genes,       file=paste(processed_data, "SP_genes.txt", sep=""),                       sep="\n", col.names=T,  row.names=F, quote=F)



