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


# load packages
library(reshape2, quietly=T, warn.conflicts=F)
library(stringr,  quietly=T, warn.conflicts=F)
library(dplyr,    quietly=T, warn.conflicts=F)
library(pipeR,    quietly=T, warn.conflicts=F)
library(psych,    quietly=T, warn.conflicts=F)
library(ggplot2,  quietly=T, warn.conflicts=F)


# directories
original_data  <- "/your/data/directory/original_data/"
processed_data <- "/your/data/directory/processed_data/"
fig            <- "/your/data/directory/fig/"
stat           <- "/your/data/directory/statistics/"
scripts        <- "/your/data/directory/scripts/"
source(paste(scripts, "plotting_parameters.R", sep=""))


# data import
design <-           read.table(paste(processed_data, "design.txt", sep=""),           sep="\t", header=T, stringsAsFactors=F)
Z      <- as.matrix(read.table(paste(processed_data, "modified_Zscores.txt", sep=""), sep="\t", header=T, stringsAsFactors=F, row.names=1))
SP     <-           read.table(paste(processed_data, "SP_genes.txt", sep=""),         sep="\t", header=T, stringsAsFactors=F)



design$genotype <- factor(design$genotype, levels=genotype$names)
idx <- order(design$genotype)
design <- design[idx,]

idx <- match(design$SampleID_, colnames(Z))
Z <- Z[, idx]

idx <- rownames(Z) %in% SP[,1]
Z <- Z[idx, ]


# boxplot
cor <- corr.test(Z)
r <- melt(cor$r)

idx <- match(r$Var1, design$SampleID_)
r$Var1_geno <- design$genotype[idx]
r$Var1_t    <- design$timepoint[idx]
r$Var1_rep  <- design$rep[idx]

idx <- match(r$Var2, design$SampleID_)
r$Var2_geno <- design$genotype[idx]
r$Var2_t    <- design$timepoint[idx]
r$Var2_rep  <- design$rep[idx]


# refine
idx <- r$Var1_rep == r$Var2_rep
r <- r[idx,]

idx <- r$Var1_t == r$Var2_t
r <- r[idx,]

idx <- r$Var1_geno != r$Var2_geno
r <- r[idx,]


# sort
r$Var1_geno <- factor(r$Var1_geno, levels=genotype$names)
r$Var2_geno <- factor(r$Var2_geno, levels=genotype$names)
levels(r$Var2_geno) <- geno.label
r$Var1_t <- factor(r$Var1_t, timepoint$names)


# plot
p <- ggplot(r, aes(x=Var1_geno, y=value, shape=Var1_t)) +
	geom_point(position=position_jitterdodge(), colour=c_black) +
	geom_boxplot(fill=c_white, outlier.shape=NA) +
	scale_shape_manual(values=timepoint$shapes) +
	facet_wrap( ~ Var2_geno, nrow=2, scales="free_x", labeller=label_parsed) +
	labs(x="", y="PCC(Z)", shape="Timepoint") +
	theme_RTN +
	theme(axis.text.x=element_text(angle=30, hjust=1, vjust=1),
		legend.position="top") +
	guides(shape=guide_legend(override.aes = list(fill=NA)))
ggsave(p, file=paste(fig, "COR_bt_genotype.boxplot.SP.pdf", sep=""), bg="transparent", width=8, height=6)


