
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
library(reshape2,   quietly=T, warn.conflicts=F)

# directories
original_data  <- "/your/data/directory/original_data/"
processed_data <- "/your/data/directory/processed_data/"
fig            <- "/your/data/directory/fig/"
stat           <- "/your/data/directory/statistics/"


# median-centered z scores function
transZ <- function(x) {
	m   <- median(x)
	mad <- mad(x, constant=1.4826)	# median absolute deviation
	meanAD <- mean(abs(x-mean(x)))		# mean absolute deviation
	if(mad == 0) {
		y <- (x-m)/(1.253314*meanAD)
	} else {
		y <- (x - m)/mad
	}
	return(y)
}


# data import
LFQ_melt      <- read.table(paste(processed_data, "LFQ_melt_table_corrected.txt", sep=""),       header=T, stringsAsFactors=F)  
design        <- read.table(paste(processed_data,  "design.txt", sep=""),                        header=T, stringsAsFactors=F)



## using raw LFQ values
# dcast
idx <- names(LFQ_melt) %in% c("SampleID", "LFQ")
dcast <- dcast(LFQ_melt[,!idx], ID ~ genotype + timepoint + rep, value.var="mLFQ")
mat <- as.matrix(dcast[,-1])
colnames(mat) <- names(dcast)[-1]
rownames(mat) <- dcast$ID

Z <- t(apply(mat, 1, transZ))

# saturate outliers
satZ <- Z

max <- quantile(satZ, .99)
idx <- satZ > max
satZ[idx] <- max

min <- quantile(satZ, .01)
idx <- satZ < min
satZ[idx] <- min

# export
write.table(Z,    file=paste(processed_data, "modified_Zscores.txt", sep=""), sep="\t", quote=F, row.names=T, col.names=NA)
write.table(satZ, file=paste(processed_data, "saturated_mod_Z.txt", sep=""),  sep="\t", quote=F, row.names=T, col.names=NA)



