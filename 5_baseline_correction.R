
#
# R script for processing data from proteome mass spec analysis
# to define lowest value
#
# orignally by Ryohei Thomas Nakano, PhD
# nakano@mpipz.mpg.de
#

options(warn=-1)

# cleaning up
rm(list=ls())

# libraries
library(ggplot2,   quietly=T, warn.conflicts=F)
library(reshape2,  quietly=T, warn.conflicts=F)
library(dplyr,     quietly=T, warn.conflicts=F)
library(pipeR,     quietly=T, warn.conflicts=F)


# directories
original_data  <- "/your/data/directory/original_data/"
processed_data <- "/your/data/directory/processed_data/"
fig            <- "/your/data/directory/fig/"
stat           <- "/your/data/directory/statistics/"
scripts        <- "/your/data/directory/scripts/"
source(paste(scripts, "plotting_parameters.R", sep=""))

# data import
LFQ_melt <- read.table(paste(processed_data, "LFQ_melt_table.txt", sep=""), header=T, stringsAsFactors=F)
design   <- read.table(paste(processed_data, "design.txt", sep=""), sep="\t", header=T, stringsAsFactors=F)




# subset to non-zero data
idx <- LFQ_melt$LFQ > 0
melt <- LFQ_melt[idx,]

# proportion of proteins that are less than particular values
idx <- order(melt$LFQ)
melt <- melt[idx,]

melt$prop <- (1:nrow(melt))/nrow(melt)

# cutoff at 1%
idx <- which.min(abs(melt$prop - 0.01))
cutoff <- melt$LFQ[idx]
cutoff_round <- round(log2(cutoff), digits=1)

# plotting
p <- ggplot(melt, aes(x=log2(LFQ), y=prop)) +
       geom_line() +
       geom_hline(yintercept=0.01, linetype="dashed", colour=c_grey) +
       geom_vline(xintercept=log2(cutoff), linetype="dotted", colour=c_cudo_magenta, size=1) +
       scale_x_continuous(breaks=c(15, 20, cutoff_round, 25, 30, 35, 40), labels=c(15,20,"",25,30,35,40)) +
       scale_y_continuous(breaks=c(0, .01, .2, .4, .6, .8, 1), labels=c("0%", "", "20%", "40%", "60%", "80%", "100%")) +
       theme_RTN +
       labs(x="log2(LFQ)", y="Proportion of proteins")
ggsave(p, file=paste(fig, "Cumulative_proportion.pdf", sep=""), bg="transparent", width=4, height=3)


# baseline correction
LFQ_melt$mLFQ <- LFQ_melt$LFQ
idx <- LFQ_melt$LFQ < cutoff
LFQ_melt$mLFQ[idx] <- cutoff

# remove proteins of which LFQs were less than the cutoff in all samples
summary <- LFQ_melt %>>% group_by(ID) %>>% summarise(count=sum(mLFQ==cutoff))

idx <- summary$count == nrow(design)
IDs <- summary$ID[idx]

idx <- LFQ_melt$ID %in% IDs
LFQ_melt <- LFQ_melt[!idx,]

# dcast
LFQ <- as.matrix(dcast(LFQ_melt[, c("ID", "genotype", "timepoint", "rep", "mLFQ")], ID ~ genotype + timepoint + rep, value.var="mLFQ"))


write.table(LFQ_melt, paste(processed_data, "LFQ_melt_table_corrected.txt", sep=""), quote=F, row.names=F, col.names=T, sep="\t")
write.table(LFQ,      paste(processed_data, "LFQ_table_corrected.txt", sep=""),      quote=F, row.names=F, col.names=T, sep="\t")
write.table(c(cutoff=cutoff, log2=log2(cutoff)), paste(stat, "baseline_correction.txt", sep=""), quote=F, row.names=T, col.names=F, sep="\t")



