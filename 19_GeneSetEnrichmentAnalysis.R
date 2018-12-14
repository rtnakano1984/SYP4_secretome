
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
library(org.At.tair.db,     quietly=T, warn.conflicts=F)
library(clusterProfiler,    quietly=T, warn.conflicts=F)
library(dplyr,              quietly=T, warn.conflicts=F)
library(ggplot2,            quietly=T, warn.conflicts=F)

# directories
original_data  <- "/your/data/directory/original_data/"
processed_data <- "/your/data/directory/processed_data/"
fig            <- "/your/data/directory/fig/"
stat           <- "/your/data/directory/statistics/"
scripts        <- "/your/data/directory/scripts/"
source(paste(scripts, "plotting_parameters.R", sep=""))


# data import
design       <-           read.table(paste(processed_data, "design.txt", sep=""),                    sep="\t", header=T,              check.names=F, stringsAsFactors=F)
logFC_P      <- as.matrix(read.table(paste(processed_data, "logFC.P.txt", sep=""),                   sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F))
SP_genes     <-           read.table(paste(processed_data, "SP_genes.txt", sep=""),                  sep="\t", header=T,              check.names=F, stringsAsFactors=F)



# subset to SP genes
idx <- rownames(logFC_P) %in% SP_genes[,1]
logFC_P <- logFC_P[idx, ]


# contrasts
contrasts <- unique(design$contrast)
contrasts <- contrasts[!is.na(contrasts)]


gsego_list <- lapply(contrasts, function(x) {

	message(paste("\n\n GSEA for ", x, "\n", sep=""))

	# getting gene list - ranked list metrics, names to be gene names (AGI codes)
	geneList <- logFC_P[, paste(x, "_logFC", sep="")]

	# decreasing sort - necessary for GSEA
	idx <- order(geneList, decreasing=T)
	geneList <- geneList[idx]

	# GSEA
	gsego <- gseGO(geneList    = geneList,
	              OrgDb        = org.At.tair.db,
	              keyType      = "TAIR",
	              ont          = "ALL",
	              nPerm        = 2500,
	              minGSSize    = 10,
	              maxGSSize    = 5000,
	              pvalueCutoff = 1,
	              verbose      = TRUE)

	summary <- summary(gsego)
	summary$contrasts <- x

	return(summary)
})

gsego_all <- rbind_all(gsego_list)

write.table(gsego_all, file=paste(stat, "GSEO_all.txt", sep=""), quote=F, sep="\t", col.names=T, row.names=F)



