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
Z            <-  as.matrix(read.table(paste(processed_data, "saturated_mod_Z.txt", sep=""),           sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F))
LFQ_melt     <-            read.table(paste(processed_data, "LFQ_melt_table_corrected.txt", sep=""),  sep="\t", header=T,              check.names=F, stringsAsFactors=F)
total        <-            read.table(paste(processed_data, "totalPeptide_melt_table.txt", sep=""),   sep="\t", header=T,              check.names=F, stringsAsFactors=F)
unique       <-            read.table(paste(processed_data, "uniquePeptide_melt_table.txt", sep=""),  sep="\t", header=T,              check.names=F, stringsAsFactors=F)
SP_genes     <-            read.table(paste(processed_data, "SP_genes.txt", sep=""),                  sep="\t", header=T,              check.names=F, stringsAsFactors=F)
sorted_logFC <-  as.matrix(read.table(paste(processed_data, "sorted_logFC.kmeans.txt", sep=""),       sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F))
k            <- as.numeric(read.table(paste(stat, "AIC_BIC_best_k.logFC.txt", sep=""),                sep="\t", header=F,              check.names=F, stringsAsFactors=F))
clusters     <-            read.table(paste(stat, "k_means_", k, ".sorted_hclust.txt", sep=""),       sep="\t", header=T,              check.names=F, stringsAsFactors=F)

sorted_ID <- rownames(sorted_logFC)


message("Processing data...")


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



# SP presense/absense table
SP <- data.frame(
	ID=factor(sorted_ID, levels=sorted_ID),
	SP=as.factor(as.numeric(sorted_ID %in% SP_genes[,1])))




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
ggsave(p, file=paste(fig, "heatmap_Z.png", sep=""), bg="transparent", width=6, height=8)


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
ggsave(p, file=paste(fig, "heatmap_LFQ.png", sep=""), bg="transparent", width=6, height=8)




# plot total peptides count
p <- ggplot(total, aes(x=SampleID, y=ID, fill=log2(Total_peptides))) +
	geom_raster() +
	scale_fill_gradient2(low="#0068B7", mid=c_white, high=c_red, na.value=c_grey) +
	labs(x="", y="", fill="Total Peptides") +
	facet_grid(. ~ genotype, space="free", scales="free", drop=T, labeller=label_parsed) +
	theme_RTN +
	theme(legend.position="top",
		strip.text=element_text(size=6),
		legend.key=element_rect(colour=c_black),
		axis.text=element_blank(),
		panel.background=element_rect(colour=c_black, size=.5, fill=NA),
		axis.line=element_blank(),
		axis.ticks=element_blank())
ggsave(p, file=paste(fig, "heatmap_total.png", sep=""), bg="transparent", width=6, height=8)



# plot unique peptides count
p <- ggplot(unique, aes(x=SampleID, y=ID, fill=log2(RazorUnique_peptides))) +
	geom_raster() +
	scale_fill_gradient2(low="#0068B7", mid=c_white, high=c_red, na.value=c_grey) +
	labs(x="", y="", fill="Razer + Unique Peptides") +
	facet_grid(. ~ genotype, space="free", scales="free", drop=T, labeller=label_parsed) +
	theme_RTN +
	theme(legend.position="top",
		strip.text=element_text(size=6),
		legend.key=element_rect(colour=c_black),
		axis.text=element_blank(),
		panel.background=element_rect(colour=c_black, size=.5, fill=NA),
		axis.line=element_blank(),
		axis.ticks=element_blank())
ggsave(p, file=paste(fig, "heatmap_unique.png", sep=""), bg="transparent", width=6, height=8)



# plot SPs
p <- ggplot(SP, aes(x=1, y=ID, fill=SP)) +
	geom_raster() +
	scale_fill_manual(values=c("0"="white", "1"=c_red)) +
	labs(x="", y="") +
	theme_RTN +
	theme(legend.position="none",
		axis.text=element_blank(),
		axis.line=element_blank(),
		axis.ticks=element_blank())
ggsave(p, file=paste(fig, "heatmap_SP.png", sep=""), bg="transparent", width=0.5, height=8)


