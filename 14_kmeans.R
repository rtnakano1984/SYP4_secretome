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
library(stringr,   quietly=T, warn.conflicts=F)
library(dplyr,     quietly=T, warn.conflicts=F)
library(pipeR,     quietly=T, warn.conflicts=F)
library(reshape2,  quietly=T, warn.conflicts=F)
library(ggplot2,   quietly=T, warn.conflicts=F)


# directories
original_data  <- "/your/data/directory/original_data/"
processed_data <- "/your/data/directory/processed_data/"
fig            <- "/your/data/directory/fig/"
stat           <- "/your/data/directory/statistics/"
scripts        <- "/your/data/directory/scripts/"
source(paste(scripts, "plotting_parameters.R", sep=""))


# function
represent <- function(fit, x) {
	ids <- names(fit$cluster[fit$cluster==x])
	
	idx <- rownames(logFC) %in% ids
	logFC_temp <- logFC[idx, ]

	if(length(ids)==1){ med <- logFC_temp } else { med <- apply(logFC_temp, 2, median) }

	return(med)
}

# data import
design  <-             read.table(paste(processed_data, "design.txt", sep=""),         sep="\t", header=T,              check.names=F, stringsAsFactors=F)
logFC   <-   as.matrix(read.table(paste(processed_data, "logFC_mat.txt", sep=""),   sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F))
k       <-  as.numeric(read.table(paste(stat, "AIC_BIC_best_k.logFC.txt", sep=""),  sep="\t", header=F,              check.names=F, stringsAsFactors=F))
 

# k means Clustering
fit <- kmeans(logFC, centers=k)

# plot individual clusters
represent_mat <- sapply(1:k, represent, fit=fit)

# corerlation between clusters
colnames(represent_mat) <- 1:k
cor <- cor(represent_mat)

# hclsut of clusters
d <- as.dist(1-cor)
hclust <- hclust(d, "average")
sorted_clusters <- rownames(cor)[hclust$order]

pdf(file=paste(fig, "dendrogram_hclust_of_k_mean_clusters_", k, ".pdf", sep=""))
plot(hclust)
dev.off()

#export sorted Z
melt <- melt(logFC)
melt$Var1 <- as.character(melt$Var1)
melt$Var2 <- str_replace(melt$Var2, "_logFC", "")
idx <- match(melt$Var2, design$contrast)
melt <- data.frame(melt, design[idx, c("genotype", "timepoint")], stringsAsFactors=F)

idx <- match(melt$Var1, names(fit$cluster))
melt$cluster <- factor(fit$cluster[idx], levels=sorted_clusters)

idx <- order(melt$cluster)
sorted_ids <- unique(melt$Var1[idx])

sorted_logFC   <- logFC[sorted_ids, ]
write.table(sorted_logFC, file=paste(processed_data, "sorted_logFC.kmeans.txt", sep=""), row.names=T, col.names=NA, sep="\t", quote=F)

# export clsuter info
cluster <- data.frame(ID=names(fit$cluster), Cluster=fit$cluster)
cluster$Cluster <- factor(cluster$Cluster, levels=sorted_clusters)
idx <- order(cluster$Cluster)
cluster <- cluster[idx,]
write.table(cluster, file=paste(stat, "k_means_", k, ".sorted_hclust.txt", sep=""), sep="\t", quote=F, col.names=T, row.names=F)

