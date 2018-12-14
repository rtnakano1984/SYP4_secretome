

#
# R script for processing data from proteome mass spec analysis
# after Signal Peptide prediction
#
# orignally by Ryohei Thomas Nakano, PhD
# nakano@mpipz.mpg.de
#

options(warn=-1)

# clean up
rm(list=ls())

# packages
library(stringr,  quietly=T, warn.conflicts=F)
library(reshape2, quietly=T, warn.conflicts=F)
library(pipeR,    quietly=T, warn.conflicts=F)
library(ggplot2,  quietly=T, warn.conflicts=F)

# directories
original_data  <- "/your/data/directory/original_data/"
processed_data <- "/your/data/directory/processed_data/"
fig            <- "/your/data/directory/fig/"
stat           <- "/your/data/directory/statistics/"
scripts        <- "/your/data/directory/scripts/"
source(paste(scripts, "plotting_parameters.R", sep=""))


# data import
design  <- read.table(paste(processed_data, "design.txt", sep=""),                      sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F)
logFC_P <- read.table(paste(processed_data, "logFC.P.txt", sep=""),                     sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F)
DE      <- read.table(paste(stat, "GLM_LRT.whole_gene_significance_table.txt", sep=""), sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F)
SP      <- read.table(paste(processed_data, "SP_genes.txt", sep=""),                    sep="\t", header=T, stringsAsFactors=F)


# extracting DEGs
idx_P  <- rowSums(abs(DE)) != 0
idx_FC <- rowSums(abs(logFC_P[, str_detect(colnames(logFC_P), "logFC")]) > log2(FC_threshold)) != 0

DEG <- intersect(rownames(DE)[idx_P], rownames(logFC_P)[idx_FC])
DEG_table <- logFC_P[match(DEG, rownames(logFC_P)), ]
write.table(DEG_table, file=paste(stat, "DEG.txt", sep=""), sep="\t", quote=F, row.names=T, col.names=NA)


# volcano plots
idx <- match(rownames(logFC_P), rownames(DE))
logFC_P <- data.frame(ID=rownames(logFC_P), logFC_P)
melt <- melt(logFC_P, id.vars="ID")
	
idx <- str_detect(melt$variable, "logFC")
melt$type[idx]  <- "logFC"
melt$type[!idx] <- "PValue"

melt$variable <- str_replace(melt$variable, "_logFC|_PValue", "")
idx <- match(melt$variable, design$contrast)
melt <- data.frame(melt, design[idx, c("genotype", "timepoint")], stringsAsFactors=F)
rownames(melt) <- NULL

dcast <- dcast(melt[, -2], ID + genotype + timepoint ~ type, value.var="value")

sig <- (dcast$PValue < alpha) & (abs(dcast$logFC) > log2(FC_threshold))
SP_idx <- dcast$ID %in% SP[,1]

dcast$sig <- 0
dcast$sig[sig & !SP_idx] <- 1
dcast$sig[sig & SP_idx] <- 2
dcast$sig <- factor(dcast$sig, levels=c(0,1,2))


dcast$timepoint <- factor(dcast$timepoint, levels=c(0,24,48))
dcast$genotype  <- factor(dcast$genotype,  levels=genotype$names)
levels(dcast$genotype) <- geno.label

p <- ggplot(dcast, aes(x=logFC, y=-log10(PValue), shape=timepoint, colour=sig)) +
	geom_point(alpha=.25) +
	geom_vline(xintercept=log2(c(1/FC_threshold, FC_threshold)), size=1, linetype="dashed", colour=c_grey) +
	geom_hline(yintercept=-log10(alpha), size=1, linetype="dashed", colour=c_grey) +
	facet_wrap( ~ genotype, scale="free", switch="x", nrow=1, labeller=label_parsed) +
	scale_shape_manual(values=c(8,1,17)) +
	scale_colour_manual(values=c(c_black, c_green, c_cudo_magenta), labels=c("n.s.", "Sig. no_SP", "Sig. SP")) +
	labs(x="log2 Fold Change", y="FDR-corrected P values (-log10)", shape="Timepoint", colour="Significance") +
	theme_RTN +
	theme(legend.position="top") +
	guides(shape=guide_legend(override.aes=list(size=3)), colour=guide_legend(override.aes=list(size=2, alpha=1)))
ggsave(p, file=paste(fig, "volcaneo_facet.png", sep=""), width=10, height=3.5, bg="transparent")



