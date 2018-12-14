

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


# directories
original_data  <- "/your/data/directory/original_data/"

# data import
dat  <- read.table(paste(original_data, "proteinGroups.txt",  sep=""), sep="\t", header=T, stringsAsFactors=F)

# extract Majority protein ID
IDs <- dat$Majority.protein.IDs
ID_list <- str_split(IDs, "\\;")

# merge splice variants
ID_list_AGI <- lapply(ID_list, function(x) str_replace(x, "\\..*", ""))

# count loci within protein groups
ID_list_count <- sapply(ID_list_AGI, function(x) length(unique(x[str_detect(x, "^AT")])))

# list of genes with progein_group_id
IDs_merged <- sapply(ID_list_AGI, function(x) paste(unique(x), collapse=";"))

out <- data.frame(ID=paste("ProteinGroup_", 1:length(IDs), sep=""), IDs=IDs, loci=IDs_merged, count=ID_list_count, stringsAsFactors=F)
write.table(out, file=paste(original_data, "Majority_protein_groups.txt", sep=""), sep="\t", quote=F, row.names=F, col.names=T)


# list of genes with multiple loci within progein groups
idx <- out$count > 1
target <- out$ID[idx]

AGI_list <- lapply(target, function(x){
			   idx <- out$ID==x
			   AGIs <- out$loci[idx]
			   split <- str_split(AGIs, ";")[[1]]
			   temp <- data.frame(ID=x, loci=split, stringsAsFactors=F)
})

target_out <- do.call(rbind, AGI_list)
write.table(target_out, file=paste(original_data, "Majority_protein_groups.with_multiple_loci.txt", sep=""),
	    sep="\t", row.names=F, col.names=T, quote=F)

idx <- str_detect(target_out$loci, "^AT") 
write.table(target_out[idx, ], file=paste(original_data, "Majority_protein_groups.with_multiple_loci.At.txt", sep=""),
	    sep="\t", row.names=F, col.names=T, quote=F)

loci_out <- target_out[idx, "loci"]
loci_out <- str_replace(tolower(loci_out), "at", "At")
write.table(loci_out, file=paste(original_data, "Majority_protein_groups.with_multiple_loci.At.loci_only.txt", sep=""),
	    sep="\t", row.names=F, col.names=F, quote=F)


