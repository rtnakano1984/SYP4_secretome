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

# packages
library(reshape2,  quietly=T, warn.conflicts=F)
library(stringr,   quietly=T, warn.conflicts=F)
library(ggplot2,   quietly=T, warn.conflicts=F)


# directories
original_data  <- "/your/data/directory/original_data/"
processed_data <- "/your/data/directory/processed_data/"
fig            <- "/your/data/directory/fig/"
stat           <- "/your/data/directory/statistics/"
scripts        <- "/your/data/directory/scripts/"
source(paste(scripts, "plotting_parameters.R", sep=""))



# data import
message("Importing data...")
design       <-            read.table(paste(processed_data, "design.txt", sep=""),                    sep="\t", header=T,              check.names=F, stringsAsFactors=F)
logFC_P      <-  as.matrix(read.table(paste(processed_data, "logFC.P.txt", sep=""),                   sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F))
Z            <-  as.matrix(read.table(paste(processed_data, "saturated_mod_Z.txt", sep=""),           sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F))
LFQ_melt     <-            read.table(paste(processed_data, "LFQ_melt_table_corrected.txt", sep=""),  sep="\t", header=T,              check.names=F, stringsAsFactors=F)
SP_genes     <-            read.table(paste(processed_data, "SP_genes.txt", sep=""),                  sep="\t", header=T,              check.names=F, stringsAsFactors=F)
total        <-            read.table(paste(processed_data, "totalPeptide_melt_table.txt", sep=""),   sep="\t", header=T,              check.names=F, stringsAsFactors=F)
unique       <-            read.table(paste(processed_data, "uniquePeptide_melt_table.txt", sep=""),  sep="\t", header=T,              check.names=F, stringsAsFactors=F)
sorted_logFC <-  as.matrix(read.table(paste(processed_data, "sorted_logFC.kmeans.SP.txt", sep=""),    sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F))
k            <- as.numeric(read.table(paste(stat, "AIC_BIC_best_k.logFC.SP.txt", sep=""),             sep="\t", header=F,              check.names=F, stringsAsFactors=F))
clusters     <-            read.table(paste(stat, "k_means_", k, ".sorted_hclust.SP.txt", sep=""),    sep="\t", header=T,              check.names=F, stringsAsFactors=F)

sorted_ID <- rownames(sorted_logFC)


message("Processing data...")
# subset to SP genes
idx <- rownames(logFC_P) %in% SP_genes[,1]
logFC_P <- logFC_P[idx, ]

idx <- rownames(Z) %in% SP_genes[,1]
Z <- Z[idx, ]

idx <- LFQ_melt$ID %in% SP_genes[,1]
LFQ_melet <- LFQ_melt[idx, ]

idx <- total$ID %in% SP_genes[,1]
total <- total[idx, ]

idx <- unique$ID %in% SP_genes[,1]
unique <- unique[idx, ]



# split logFC_P
idx <- str_detect(colnames(logFC_P), "logFC")
logFC <- logFC_P[, idx]
P     <- logFC_P[, !idx]

# saturate logFC
max <- quantile(logFC, .99)
idx <- (logFC > max)
logFC[idx] <- max

min <- quantile(logFC, .01)
idx <- (logFC < min)
logFC[idx] <- min



# melt and sort logFC
melt_FC <- melt(logFC)
melt_FC$Var1 <- factor(melt_FC$Var1, levels=sorted_ID)
melt_FC$Var2 <- str_replace(melt_FC$Var2, "_logFC", "")

idx <- match(melt_FC$Var2, design$contrast)
melt_FC <- data.frame(melt_FC, design[idx,], stringsAsFactors=F)

melt_FC$SampleID  <- factor(melt_FC$SampleID,  levels=design$SampleID)
melt_FC$timepoint <- factor(melt_FC$timepoint, levels=timepoint$names)
melt_FC$genotype  <- factor(melt_FC$genotype,  levels=genotype$names)
levels(melt_FC$genotype) <- geno.label


# melt and sort P
melt_P <- melt(P)
melt_P$Var1 <- factor(melt_P$Var1, levels=sorted_ID)
melt_P$Var2 <- str_replace(melt_P$Var2, "_PValue", "")

idx <- match(melt_P$Var2, design$contrast)
melt_P <- data.frame(melt_P, design[idx,], stringsAsFactors=F)

melt_P$SampleID  <- factor(melt_P$SampleID,  levels=design$SampleID)
melt_P$timepoint <- factor(melt_P$timepoint, levels=timepoint$names)
melt_P$genotype  <- factor(melt_P$genotype,  levels=genotype$names)
levels(melt_P$genotype) <- geno.label

idx_P     <- melt_P$value < alpha
sigID_P <- paste(melt_P$Var1, melt_P$Var2, sep="_")[idx_P]

idx_logFC <- melt_FC$value < FC_threshold
sigID_logFC <- paste(melt_FC$Var1, melt_FC$Var2, sep="_")[idx_logFC]

sigID <- intersect(sigID_P, sigID_logFC)
idx <- paste(melt_P$Var1, melt_P$Var2, sep="_") %in% sigID
melt_P$sig <- as.numeric(idx)


# melt and sort Z
melt_Z <- melt(Z)
melt_Z$Var1 <- factor(melt_Z$Var1, levels=sorted_ID)
melt_Z$Var2 <- str_replace(melt_Z$Var2, "_logFC", "")

idx <- match(melt_Z$Var2, design$SampleID_)
melt_Z <- data.frame(melt_Z, design[idx,], stringsAsFactors=F)

melt_Z$SampleID_  <- factor(melt_Z$SampleID_, levels=design$SampleID_)
melt_Z$timepoint  <- factor(melt_Z$timepoint, levels=timepoint$names)
melt_Z$genotype   <- factor(melt_Z$genotype,  levels=genotype$names)
levels(melt_Z$genotype) <- geno.label



# sort LFQ_melt
LFQ_melt$SampleID_  <- paste(LFQ_melt$genotype, LFQ_melt$timepoint, LFQ_melt$rep, sep="_")
idx <- match(LFQ_melt$SampleID_, design$SampleID_)
LFQ_melt[, c("genotype", "timepoint", "rep")] <- design[idx, c("genotype", "timepoint", "rep")]

LFQ_melt$ID         <- factor(LFQ_melt$ID,        levels=sorted_ID)
LFQ_melt$SampleID_  <- factor(LFQ_melt$SampleID_, levels=design$SampleID_)
LFQ_melt$timepoint  <- factor(LFQ_melt$timepoint, levels=timepoint$names)
LFQ_melt$genotype   <- factor(LFQ_melt$genotype,  levels=genotype$names)
levels(LFQ_melt$genotype) <- geno.label


# sort clusters
clusters$ID <- factor(clusters$ID, levels=sorted_ID)

idx <- order(clusters$ID)
clusters <- clusters[idx,]

clusters$Cluster <- factor(clusters$Cluster, levels=unique(clusters$Cluster))
clusters$Cluster_bin <- as.numeric(clusters$Cluster)
clusters$Cluster_bin <- as.factor(clusters$Cluster_bin %% 2)



# sort total
total$SampleID <- str_replace_all(total$SampleID, "Peptides.", "")
total$genotype <- str_replace(total$genotype, "Col", "Col-0")
total$genotype <- str_replace(total$genotype, "pen1", "pen1-1")

total$timepoint <- str_replace(total$timepoint, "h", "")

total$SampleID   <- factor(total$SampleID, levels=design$SampleID)
total$ID         <- factor(total$ID, levels=sorted_ID)
total$timepoint  <- factor(total$timepoint, levels=timepoint$names)
total$genotype   <- factor(total$genotype,  levels=genotype$names)
levels(total$genotype) <- geno.label


# sort unique
unique$SampleID <- str_replace_all(unique$SampleID, "Razor...unique.peptides.", "")
unique$genotype <- str_replace(unique$genotype, "Col", "Col-0")
unique$genotype <- str_replace(unique$genotype, "pen1", "pen1-1")

unique$timepoint <- str_replace(unique$timepoint, "h", "")

unique$SampleID   <- factor(unique$SampleID, levels=design$SampleID)
unique$ID         <- factor(unique$ID, levels=sorted_ID)
unique$timepoint  <- factor(unique$timepoint, levels=timepoint$names)
unique$genotype   <- factor(unique$genotype,  levels=genotype$names)
levels(unique$genotype) <- geno.label




message("Plotting data...")

# plot logFC
p <- ggplot(melt_FC, aes(x=Var2, y=Var1, fill=value)) +
	geom_raster() +
	scale_fill_gradient2(high=c_dark_green, mid="white", low=c_cudo_magenta, midpoint=0, na.value="gray95") +
	labs(x="", y="", fill="logFC") +
	facet_grid(. ~ genotype, space="free", scales="free", drop=T, labeller=label_parsed) +
	theme_RTN +
	theme(legend.position="top",
		strip.text=element_text(size=6),
		axis.text=element_blank(),
		axis.line=element_blank(),
		panel.background=element_rect(colour=c_black, size=.5, fill=NA),
		axis.ticks=element_blank()) +
	guides(fill=guide_colourbar(barheight=1, barwidth=5))
ggsave(p, file=paste(fig, "heatmap_logFC.SP.png", sep=""), bg="transparent", width=6, height=8)

p <- p + scale_fill_gradient2(low=c_blue, mid=c_white, high=c_red, midpoint=0, na.value="grey95")
ggsave(p, file=paste(fig, "heatmap_logFC.SP.blue_red.png", sep=""), bg="transparent", width=6, height=8)


# plot P
p <- ggplot(melt_P, aes(x=Var2, y=Var1, fill=factor(sig))) +
	geom_raster() +
	scale_fill_manual(values=c("0"="white", "1"=c_red)) +
	labs(x="", y="", fill="Significance") +
	facet_grid(. ~ genotype, space="free", scales="free", drop=T, labeller=label_parsed) +
	theme_RTN +
	theme(legend.position="top",
		strip.text=element_text(size=6),
		legend.key=element_rect(colour=c_black),
		axis.text=element_blank(),
		panel.background=element_rect(colour=c_black, size=.5, fill=NA),
		axis.line=element_blank(),
		axis.ticks=element_blank())
ggsave(p, file=paste(fig, "heatmap_P.SP.png", sep=""), bg="transparent", width=6, height=8)



# plot whole Z
p <- ggplot(melt_Z, aes(x=Var2, y=Var1, fill=value)) +
	geom_raster() +
	scale_fill_gradient2(low="#0068B7", mid="white", high="#F39800", midpoint=0) +
	labs(x="", y="", fill="Modified Z Score") +
	facet_grid(. ~ genotype, space="free", scales="free", drop=T, labeller=label_parsed) +
	theme_RTN +
	theme(legend.position="top",
		strip.text=element_text(size=6),
		legend.key=element_rect(colour=c_black),
		axis.text=element_blank(),
		panel.background=element_rect(colour=c_black, size=.5, fill=NA),
		axis.line=element_blank(),
		axis.ticks=element_blank())
ggsave(p, file=paste(fig, "heatmap_Z.SP.png", sep=""), bg="transparent", width=6, height=8)


# plot whole LFQ
mid <- median(log2(LFQ_melt$mLFQ))
p <- ggplot(LFQ_melt, aes(x=SampleID_, y=ID, fill=log2(mLFQ))) +
	geom_raster() +
	scale_fill_gradient2(low=c_cudo_skyblue, mid=c_black, high=c_yellow, midpoint=mid) +
	labs(x="", y="", fill="LFQ") +
	facet_grid(. ~ genotype, space="free", scales="free", drop=T, labeller=label_parsed) +
	theme_RTN +
	theme(legend.position="top",
		strip.text=element_text(size=6),
		legend.key=element_rect(colour=c_black),
		axis.text=element_blank(),
		panel.background=element_rect(colour=c_black, size=.5, fill=NA),
		axis.line=element_blank(),
		axis.ticks=element_blank())
ggsave(p, file=paste(fig, "heatmap_LFQ.SP.png", sep=""), bg="transparent", width=6, height=8)


# plot clusters
p <- ggplot(clusters, aes(x=1, y=ID, fill=Cluster_bin)) +
	geom_raster() +
	scale_fill_manual(values=c(c_black, c_grey), guide=F) +
	labs(x="", y="") +
	theme_RTN +
	theme(legend.position="none",
		axis.text=element_blank(),
		axis.line=element_blank(),
		axis.ticks=element_blank())
ggsave(p, file=paste(fig, "heatmap_cluster.SP.png", sep=""), bg="transparent", width=0.5, height=8)



