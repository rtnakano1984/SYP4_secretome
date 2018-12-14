
#
# R script for processing data from proteome mass spec analysis
#
# orignally by Ryohei Thomas Nakano, PhD
# nakano@mpipz.mpg.de
#

options(warn=-1)

# cleaning up

rm(list=ls())


# loading libraries

library(stringr,   quietly=T, warn.conflicts=F)
library(reshape2,  quietly=T, warn.conflicts=F)


# directories
original_data  <- "/your/data/directory/original_data/"
processed_data <- "/your/data/directory/processed_data/"
fig            <- "/your/data/directory/fig/"
stat           <- "/your/data/directory/statistics/"

# data import
dat  <- read.table(paste(processed_data, "proteinGroups_curated.txt",  sep=""), sep="\t", header=T, stringsAsFactors=F)
design_table  <- read.table(paste(original_data, "design.txt",  sep=""), sep="\t", header=T, stringsAsFactors=F)

# select At proteins
idx <- str_detect(dat$ID, "^AT[1-5]G")
dat <- dat[idx,]


# subset to rep1-4
idx <- design_table$rep %in% 1:4
design_table <- design_table[idx, ]

idx <- str_detect(colnames(dat), "h\\.5$")
dat <- dat[, !idx]


# subset to LFQs
idx <- str_detect(names(dat), "^LFQ")
LFQ <- data.frame(ID=dat$ID, dat[, idx])

# subset to razor + unique peptide counts
idx <- str_detect(names(dat), "^Razor...unique.*.[1-5]")
unique <- data.frame(ID=dat$ID, dat[,idx])

# subset to total peptide counts
idx <- str_detect(names(dat), "^Peptide.*.[1-5]")
total <- data.frame(ID=dat$ID, dat[,idx])

# subset to iBAQs
idx <- str_detect(names(dat), "^iBAQ.*.[1-5]")
iBAQ <- data.frame(ID=dat$ID, dat[,idx])



# melt LFQ
LFQ_melt <- melt(LFQ, id.vars="ID")
names(LFQ_melt) <- c("ID", "SampleID", "LFQ")

# replace NA with zero
idx <- is.na(LFQ_melt$LFQ)
LFQ_melt$LFQ[idx] <- 0

# assign meta info
design <- matrix(unlist(str_split(LFQ_melt$SampleID, "\\.")), ncol=5, byrow=T)
LFQ_melt$genotype  <- design[,3]
LFQ_melt$timepoint <- design[,4]
LFQ_melt$rep       <- design[,5]

# remove rep#5
idx <- LFQ_melt$rep != 5
LFQ_melt <- LFQ_melt[idx,]

# sort
idx <- order(LFQ_melt$ID, LFQ_melt$genotype, LFQ_melt$timepoint, LFQ_melt$rep)
LFQ_melt <- LFQ_melt[idx, ]

# format
LFQ_melt$LFQ  <- format(LFQ_melt$LFQ,  scientific=F)



# melt unique
unique_melt <- melt(unique, id.vars="ID")
names(unique_melt) <- c("ID", "SampleID", "RazorUnique_peptides")

# assign meta info
design <- matrix(unlist(str_split(unique_melt$SampleID, "\\.")), ncol=8, byrow=T)
unique_melt$genotype  <- design[,6]
unique_melt$timepoint <- design[,7]
unique_melt$rep       <- design[,8]

# remove rep#5
idx <- unique_melt$rep != 5
unique_melt <- unique_melt[idx,]

# sort
idx <- order(unique_melt$ID, unique_melt$genotype, unique_melt$timepoint, unique_melt$rep)
unique_melt <- unique_melt[idx, ]



# melt unique
total_melt <- melt(total, id.vars="ID")
names(total_melt) <- c("ID", "SampleID", "Total_peptides")

# assign meta info
design <- matrix(unlist(str_split(total_melt$SampleID, "\\.")), ncol=4, byrow=T)
total_melt$genotype  <- design[,2]
total_melt$timepoint <- design[,3]
total_melt$rep       <- design[,4]

# remove rep#5
idx <- total_melt$rep != 5
total_melt <- total_melt[idx,]

# sort
idx <- order(total_melt$ID, total_melt$genotype, total_melt$timepoint, total_melt$rep)
total_melt <- total_melt[idx, ]





# melt unique
iBAQ_melt <- melt(iBAQ, id.vars="ID")
names(iBAQ_melt) <- c("ID", "SampleID", "iBAQ")

# assign meta info
design <- matrix(unlist(str_split(iBAQ_melt$SampleID, "\\.")), ncol=4, byrow=T)
iBAQ_melt$genotype  <- design[,2]
iBAQ_melt$timepoint <- design[,3]
iBAQ_melt$rep       <- design[,4]

# remove rep#5
idx <- iBAQ_melt$rep != 5
iBAQ_melt <- iBAQ_melt[idx,]

# sort
idx <- order(iBAQ_melt$ID, iBAQ_melt$genotype, iBAQ_melt$timepoint, iBAQ_melt$rep)
iBAQ_melt <- iBAQ_melt[idx, ]





# output
write.table(LFQ,         file=paste(processed_data, "LFQ_table.txt", sep=""),                quote=F, col.names=T, row.names=F, sep="\t")
write.table(LFQ_melt,    file=paste(processed_data, "LFQ_melt_table.txt", sep=""),           quote=F, col.names=T, row.names=F, sep="\t")
write.table(iBAQ,        file=paste(processed_data, "iBAQ_table.txt", sep=""),               quote=F, col.names=T, row.names=F, sep="\t")
write.table(iBAQ_melt,   file=paste(processed_data, "iBAQ_melt_table.txt", sep=""),          quote=F, col.names=T, row.names=F, sep="\t")
write.table(unique_melt, file=paste(processed_data, "uniquePeptide_melt_table.txt", sep=""), quote=F, col.names=T, row.names=F, sep="\t")
write.table(total_melt,  file=paste(processed_data, "totalPeptide_melt_table.txt", sep=""),  quote=F, col.names=T, row.names=F, sep="\t")
write.table(dat$ID,      file=paste(processed_data, "ID_list.txt", sep=""),                  quote=F, col.names=F, row.names=F, sep="\n")
write.table(design_table,file=paste(processed_data, "design.txt", sep=""),                   quote=F, col.names=T, row.names=F, sep="\t")


