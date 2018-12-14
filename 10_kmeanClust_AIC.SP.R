
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
library(ggplot2,   quietly=T, warn.conflicts=F)
library(reshape2,  quietly=T, warn.conflicts=F)
library(stringr,   quietly=T, warn.conflicts=F)

# directories
original_data  <- "/your/data/directory/original_data/"
processed_data <- "/your/data/directory/processed_data/"
fig            <- "/your/data/directory/fig/"
stat           <- "/your/data/directory/statistics/"
scripts        <- "/your/data/directory/scripts/"
source(paste(scripts, "plotting_parameters.R", sep=""))


# functions
kmeansAIC <- function(fit){
	m <- ncol(fit$centers)
	n <- length(fit$cluster)
	k <- nrow(fit$centers)
	D <- fit$tot.withinss
	return(c(
		AIC=(D + 2*m*k),
		BIC=(D + log(n)*m*k)))
}

# data import
design   <-           read.table(paste(processed_data,  "design.txt",  sep=""), sep="\t", header=T,              check.names=F, stringsAsFactors=F)
logFC_P  <- as.matrix(read.table(paste(processed_data, "logFC.P.txt", sep=""),  sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F))
SP_genes <-           read.table(paste(processed_data, "SP_genes.txt", sep=""), sep="\t", header=T,              check.names=F, stringsAsFactors=F)

# subset to SP genes
idx <- rownames(logFC_P) %in% SP_genes[,1]
logFC_P <- logFC_P[idx, ]

# logFC matrix
idx <- str_detect(colnames(logFC_P), "logFC")
logFC <- logFC_P[, idx]

write.table(logFC, paste(processed_data, "logFC_mat.SP.txt", sep=""), quote=F, sep="\t", row.names=T, col.names=NA)

# k means clustering with a range of cluster numbers and calculate AIC and BIC
min <- 1
max <- 200
range <- c(min:max)

AIC <- t(sapply(range, function(x) kmeansAIC(kmeans(logFC, centers=x)) ))
AIC <- data.frame(n=range, AIC)

# plot AIC/BIC
melt <- melt(AIC, id.vars="n")
p <- ggplot(melt, aes(x=n, y=value, colour=variable)) +
	geom_line() +
	geom_point() +
	labs(colour="", x="Number of clusters", y="") +
	theme_RTN +
	theme(legend.position="top")
ggsave(p, file=paste(fig, "AIC_BIC.SP.pdf", sep=""), width=4.5, height=3.5)

# out put best number of k
k <- range[which.min(AIC$BIC)]
out <- data.frame("#AIC/BIC Likelihodd test",
	paste("#Tested range of clusters:  from ", min, " to ", max, sep=""),
	"#Best k (with lowest BIC) is:      ",
	k)
write.table(out, file=paste(stat, "AIC_BIC_best_k.logFC.SP.txt", sep=""), quote=F, row.names=F, col.names=F, sep="\n")




