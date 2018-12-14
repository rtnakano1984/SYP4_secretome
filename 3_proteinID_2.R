

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
processed_data <- "/your/data/directory/processed_data/"
fig            <- "/your/data/directory/fig/"
stat           <- "/your/data/directory/statistics/"

# data import
dat     <- read.table(paste(original_data, "proteinGroups.txt",  sep=""), sep="\t", header=T, stringsAsFactors=F)
curated <- read.table(paste(original_data, "duplicated_atted.manually-curated.txt",  sep=""), sep="\t", header=T, stringsAsFactors=F, quote="")
map     <- read.table(paste(original_data, "Majority_protein_groups.txt", sep=""), sep="\t", header=T, stringsAsFactors=F)


# merge
idx <- match(curated$ID, map$ID)
map$loci[idx] <- curated$loci

idx <- match(map$IDs, dat$Majority.protein.IDs)
dat$ID <- map$loci[idx]
dat$groupID <- map$ID[idx]


# check duplication
idx <- sapply(1:nrow(dat), function(x) {
	dat$ID[x] %in% dat$ID[-x]
	})
dup <- unique(dat$ID[idx])

for(x in dup[!is.na(dup)]){
	idx <- dat$ID == x
	temp <- dat[idx,]
	dat <- dat[!idx,]
	razor <- temp$Razor...unique.peptides
	idx <- which.max(razor)
	repres <- temp[idx,]
	dat <- rbind(dat, repres)
}

# export
write.table(dat, paste(processed_data, "proteinGroups_curated.txt", sep=""), quote=F, col.names=T, row.names=F, sep="\t")




