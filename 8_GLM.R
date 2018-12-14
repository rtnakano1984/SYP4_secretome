
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
library(edgeR,    quietly=T, warn.conflicts=F)
library(stringr,  quietly=T, warn.conflicts=F)

# directories
original_data  <- "/your/data/directory/original_data/"
processed_data <- "/your/data/directory/processed_data/"
fig            <- "/your/data/directory/fig/"
stat           <- "/your/data/directory/statistics/"
scripts        <- "/your/data/directory/scripts/"
source(paste(scripts, "plotting_parameters.R", sep=""))


# data import
message("data import")
dat_tab <- read.table(paste(processed_data, "LFQ_table_corrected.txt", sep=""), sep="\t", stringsAsFactors=F, row.names=1, header=T, check.names=F)
design  <- read.table(paste(processed_data,  "design.txt", sep=""),  sep="\t", stringsAsFactors=F,              header=T, check.names=F)

# sort data table accorindg to AGI codes
idx <- order(rownames(dat_tab))
dat_tab <- dat_tab[idx,]

idx <- match(colnames(dat_tab), design$SampleID_)
dat_tab <- dat_tab[,idx]

# DGEList object
x <- as.matrix(dat_tab)
y <- DGEList(counts=x)


# Create model
design$group <- paste(design$genotype, design$timepoint, sep="_")

rep       <- factor(design$rep,   levels=unique(design$rep))
group     <- factor(design$group, levels=unique(design$group))

model <- model.matrix( ~ 0 + group + rep)
colnames(model) <- str_replace(colnames(model), "group", "")
colnames(model) <- str_replace(colnames(model), "-", "")

message("Estimating dispersion")
y <- estimateGLMCommonDisp(y, model)
y <- estimateGLMTrendedDisp(y, model)
y <- estimateGLMTagwiseDisp(y, model)

message("GLM fitting")
fit <- glmFit(y, model)


# LRT
contrasts <- makeContrasts(
			   Col_24            = (              Col0_24 -  Col0_0   ),
			   Col_48            = (              Col0_48 -  Col0_0   ),
			   syp4243_0         = (        syp42syp43_0  -  Col0_0   ),
			   syp4243_24        = (        syp42syp43_24 -  Col0_24  ),
			   syp4243_48        = (        syp42syp43_48 -  Col0_48  ),
			   syp4243vamp721_0  = ( syp42syp43vamp721_0  -  Col0_0   ),
			   syp4243vamp721_24 = ( syp42syp43vamp721_24 -  Col0_24  ),
			   syp4243vamp721_48 = ( syp42syp43vamp721_48 -  Col0_48  ),
			   syp4243vamp722_0  = ( syp42syp43vamp722_0  -  Col0_0   ),
			   syp4243vamp722_24 = ( syp42syp43vamp722_24 -  Col0_24  ),
			   syp4243vamp722_48 = ( syp42syp43vamp722_48 -  Col0_48  ),
			   syp4243sid2_0     = (    syp42syp43sid2_0  -  Col0_0   ),
			   syp4243sid2_24    = (    syp42syp43sid2_24 -  Col0_24  ),
			   syp4243sid2_48    = (    syp42syp43sid2_48 -  Col0_48  ),
			   pen1_0            = (             pen12_0  -  Col0_0   ),
			   pen1_24           = (             pen12_24 -  Col0_24  ),
			   pen1_48           = (             pen12_48 -  Col0_48  ),
			   levels=model)
contrast_names <- colnames(contrasts)
n <- length(contrast_names)



# DEG selection
message("LRT for DEG identification")
# LRT for each contrasts
LRT.list <- lapply(1:n, function(x) glmLRT(fit, contrast=contrasts[,x]))
names(LRT.list) <- contrast_names

# logFC and PValue tables
logFC_P.list <- lapply(1:n, function(x) {
	table <- LRT.list[[x]]$table[,c(1,4)]
	table$PValue <- p.adjust(table$PValue, method=p.adj.method)
	colnames(table) <- paste(contrast_names[x], colnames(table), sep="_")
	return(table)
	})
logFC_P <- do.call(data.frame, logFC_P.list)
write.table(logFC_P, file=paste(processed_data, "logFC.P.txt", sep=""), sep="\t", quote=F, col.names=NA, row.names=T)

# Significance picking for each tested model
DE.list <- lapply(1:n, function(x) decideTestsDGE(LRT.list[[x]], adjust.method=p.adj.method, p.value=alpha))
names(DE.list) <- contrast_names

# Number of significant differentially abundant OTUs
total    <- sapply(DE.list, function(x) sum(abs(x)))
induced  <- sapply(DE.list, function(x) sum(x ==  1))
reduced  <- sapply(DE.list, function(x) sum(x == -1))
count <- data.frame(total, induced, reduced)
rownames(count) <- contrast_names
write.table(count, file=paste(stat, "number_of_DEGs.txt", sep=""), quote=F, row.names=T, col.names=NA, sep="\t")

# significance table
DE <- sapply(1:n, function(x) DE.list[[x]][,1])
colnames(DE) <- contrast_names
write.table(DE, file=paste(stat, "GLM_LRT.whole_gene_significance_table.txt", sep=""), sep="\t", quote=T, row.names=T, col.names=NA)











