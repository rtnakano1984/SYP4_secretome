
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
library(stringr,            quietly=T, warn.conflicts=F)
library(ggplot2,            quietly=T, warn.conflicts=F)

# directories
original_data  <- "/your/data/directory/original_data/"
processed_data <- "/your/data/directory/processed_data/"
fig            <- "/your/data/directory/fig/"
stat           <- "/your/data/directory/statistics/"
scripts        <- "/your/data/directory/scripts/"
source(paste(scripts, "plotting_parameters.R", sep=""))


# data import
design       <-            read.table(paste(processed_data, "design.txt", sep=""),                 sep="\t", header=T, check.names=F, stringsAsFactors=F)
k            <- as.numeric(read.table(paste(stat, "AIC_BIC_best_k.logFC.SP.txt", sep=""),          sep="\t", header=F, check.names=F, stringsAsFactors=F))
clusters     <-            read.table(paste(stat, "k_means_", k, ".sorted_hclust.SP.txt", sep=""), sep="\t", header=T, check.names=F, stringsAsFactors=F)
SP_genes     <-            read.table(paste(processed_data, "SP_genes.txt", sep=""),               sep="\t", header=T, check.names=F, stringsAsFactors=F)[,1]


# prepare cluster list
sorted_cluster <- unique(clusters$Cluster)
cluster_list <- lapply(sorted_cluster, function(x) clusters$ID[clusters$Cluster == x] )
names(cluster_list) <- paste("Cluster", sorted_cluster, sep="_")


# compare clusters vs whole genome
cg <- compareCluster(
				geneCluster=cluster_list,
				fun="enrichGO",
                OrgDb         = org.At.tair.db,
                keyType       = "TAIR",
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
                minGSSize     = 1,
                maxGSSize     = 2500)
cg_df <- as.data.frame(cg)

idx <- cg_df$p.adjust < 0.05
cg_df <- cg_df[idx,]
write.table(cg_df, file=paste(stat, "Clusterwise_enrichGO_results.whole_genome.significant_GOs.txt", sep=""), quote=F, col.names=T, row.names=F, sep="\t")



