

#
# R script for processing data from proteome mass spec analysis
#
# orignally by Ryohei Thomas Nakano, PhD
# nakano@mpipz.mpg.de
#

options(warn=-1)

# cleaning up
rm(list=ls())


# directories
original_data  <- "/your/data/directory/original_data/"

# data import
df    <- read.table(paste(original_data, "Majority_protein_groups.with_multiple_loci.txt",  sep=""), sep="\t", header=T, stringsAsFactors=F)
atted <- read.table(paste(original_data, "Majority_protein_groups.with_multiple_loci.atted.csv", sep=""), header=T, stringsAsFactors=F, sep=",")

# genes not listed in atted table
atted <- unique(atted)

idx <- df$loci %in% toupper(atted$Locus)
not_listed <- unique(df[!idx,]$ID)

idx <- df$ID %in% not_listed
out <- df[idx,]
write.table(out, file=paste(original_data, "duplicated-not_listed_atted.txt", sep=""), quote=F, sep="\t", row.names=F, col.names=T)

df <- df[!idx,]
idx <- match(df$loci, toupper(atted$Locus))
df <- data.frame(df, atted[idx,], stringsAsFactors=F)

write.table(df, file=paste(original_data, "duplicated_atted.txt", sep=""), quote=F, sep="\t", row.names=F, col.names=T)




