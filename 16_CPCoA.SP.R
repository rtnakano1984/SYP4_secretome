#
# R script for processing data from proteome mass spec analysis
# after Signal Peptide prediction
#
# orignally by Ruben Garrido-Oter, MPIPZ
# 
# editied by Ryohei Thomas Nakano, MPIPZ
# nakano@mpipz.mpg.de
#


options(warn=-1)

# clean up
rm(list=ls())


# load packages
library(ggplot2, quietly=T, warn.conflicts=F)
library(scales,  quietly=T, warn.conflicts=F)
library(grid,    quietly=T, warn.conflicts=F)
library(vegan,   quietly=T, warn.conflicts=F)


# directories
original_data  <- "/your/data/directory/original_data/"
processed_data <- "/your/data/directory/processed_data/"
fig            <- "/your/data/directory/fig/PCoA/"         #!
stat           <- "/your/data/directory/statistics/"
scripts        <- "/your/data/directory/scripts/"
source(paste(scripts, "plotting_parameters.R", sep=""))


# data import
design <- read.table(paste(processed_data, "design.txt", sep=""),  sep="\t", header=T, stringsAsFactors=F)
Z      <- read.table(paste(processed_data, "modified_Zscores.txt", sep=""), sep="\t", header=T, stringsAsFactors=F, row.names=1)
SP     <- read.table(paste(processed_data, "SP_genes.txt", sep=""), sep="\t", header=T, stringsAsFactors=F)

design <- data.frame(sapply(design, as.character))


# functions
variability_table <- function(cca){

        chi <- c(cca$tot.chi,
                       cca$CCA$tot.chi, cca$CA$tot.chi)
        variability_table <- cbind(chi, chi/chi[1])
        colnames(variability_table) <- c("inertia", "proportion")
        rownames(variability_table) <- c("total", "constrained", "unconstrained")
        return(variability_table)

}


# subset to SP proteins
idx <- rownames(Z) %in% SP[,1]
Z <- Z[idx, ]


# dissimilarity matrix
dist <- 1-cor(Z)




# normal CPoA

k <- 2
pcoa <- cmdscale(dist, k=k, eig=T)
points <- pcoa$points
eig <- pcoa$eig
points <- as.data.frame(points)
colnames(points) <- c("x", "y")

points <- cbind(points, design[match(rownames(points), design$SampleID_), ])

# assign genotypes to the corresponding colors
points$genotype    <- factor(points$genotype,    levels=genotype$names)
points$timepoint   <- factor(points$timepoint,   levels=timepoint$names)


# plot CPCo 1 and 2
p <- ggplot(points, aes(x=x, y=y, colour=genotype, shape=timepoint)) +
	stat_ellipse(type="norm", linetype=2, aes(group=genotype), alpha=.3, level=.75) +
	geom_point(alpha=.7, size=2.5) +
	scale_colour_manual(values=genotype$colours, labels=geno.label) +
	scale_shape_manual(values=timepoint$shapes) +
	labs(x=paste("PC1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
	     y=paste("PC2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
	     colour="Genotype", shape="Timepoint") +
	theme_RTN +
	theme(legend.position="top",
	      legend.text.align=0)
ggsave(p, file=paste(fig, "PCoA_corZ.whole.SP.pdf", sep=""), bg="transparent", width=7, height=4)






# CPCOA contrained by genotype

sqrt_transform <- T
d <- design

capscale.gen <- capscale(dist ~ genotype + Condition(timepoint + rep), data=d, add=F, sqrt.dist=sqrt_transform)

# ANOVA-like permutation analysis
perm_anova.gen <- anova.cca(capscale.gen)
print(perm_anova.gen)
					    
# generate variability tables and calculate confidence intervals for the variance
var_tbl.gen <- variability_table(capscale.gen)
eig <- capscale.gen$CCA$eig
variance <- var_tbl.gen["constrained", "proportion"]
p.val <- perm_anova.gen[1, 4]

# extract the weighted average (sample) scores

points <- capscale.gen$CCA$wa[, 1:2]
points <- as.data.frame(points)
colnames(points) <- c("x", "y")

points <- cbind(points, design[match(rownames(points), design$SampleID_), ])

# assign genotypes to the corresponding colors
points$genotype    <- factor(points$genotype,    levels=genotype$names)
points$timepoint   <- factor(points$timepoint,   levels=timepoint$names)

# plot CPCo 1 and 2
p <- ggplot(points, aes(x=x, y=y, colour=genotype, shape=timepoint)) +
	stat_ellipse(type="norm", linetype=2, aes(group=genotype), alpha=.3, level=.75) +
	geom_point(alpha=.7, size=2.5) +
	scale_colour_manual(values=genotype$colours, labels=geno.label) +
	scale_shape_manual(values=timepoint$shapes) +
	labs(x=paste("CPCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
	     y=paste("CPCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
	     colour="Genotype", shape="Timepoint") +
	ggtitle(paste(format(100 * variance, digits=3), " % of variance; p=",
		      format(p.val, digits=2),
		      sep="")) +
	theme_RTN +
	theme(legend.position="top",
	      legend.text.align=0)
ggsave(p, file=paste(fig, "CPCoA_genotypes.SP.pdf", sep=""), bg="transparent", width=8, height=6.5)





# CPCOA contrained by timpoint

sqrt_transform <- T
d <- design

capscale.gen <- capscale(dist ~ timepoint + Condition(genotype + rep), data=d, add=F, sqrt.dist=sqrt_transform)

# ANOVA-like permutation analysis
perm_anova.gen <- anova.cca(capscale.gen)
print(perm_anova.gen)
					    
# generate variability tables and calculate confidence intervals for the variance
var_tbl.gen <- variability_table(capscale.gen)
eig <- capscale.gen$CCA$eig
variance <- var_tbl.gen["constrained", "proportion"]
p.val <- perm_anova.gen[1, 4]

# extract the weighted average (sample) scores

points <- capscale.gen$CCA$wa[, 1:2]
points <- as.data.frame(points)
colnames(points) <- c("x", "y")

points <- cbind(points, design[match(rownames(points), design$SampleID_), ])

# assign genotypes to the corresponding colors
points$genotype    <- factor(points$genotype,    levels=genotype$names)
points$timepoint   <- factor(points$timepoint,   levels=timepoint$names)

# plot CPCo 1 and 2
p <- ggplot(points, aes(x=x, y=y, colour=genotype, shape=timepoint)) +
	stat_ellipse(type="norm", linetype=2, aes(group=timepoint), alpha=.3, level=.75) +
	geom_point(alpha=.7, size=2.5) +
	scale_colour_manual(values=genotype$colours, labels=geno.label) +
	scale_shape_manual(values=timepoint$shapes) +
	labs(x=paste("CPCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
	     y=paste("CPCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
	     colour="Genotype", shape="Timepoint") +
	ggtitle(paste(format(100 * variance, digits=3), " % of variance; p=",
		      format(p.val, digits=2),
		      sep="")) +
	theme_RTN +
	theme(legend.position="top",
	      legend.text.align=0)
ggsave(p, file=paste(fig, "CPCoA_timepoint.SP.pdf", sep=""), bg="transparent", width=8, height=6.5)





# normal PCoA and PCoA constrained by genotype for each timepoint

t <- unique(design$timepoint)
for(ti in t) {

	idx <- design$timepoint == ti
	design_temp <- design[idx,]

	idx <- colnames(Z) %in% design_temp$SampleID_
	Z_temp <- Z[, idx]

	dist <- 1-cor(Z_temp)


	sqrt_transform <- T
	d <- design_temp

	capscale.gen <- capscale(dist ~ genotype + Condition(timepoint + rep), data=d, add=F, sqrt.dist=sqrt_transform)

	# ANOVA-like permutation analysis
	perm_anova.gen <- anova.cca(capscale.gen)
	print(perm_anova.gen)
						    
	# generate variability tables and calculate confidence intervals for the variance
	var_tbl.gen <- variability_table(capscale.gen)
	eig <- capscale.gen$CCA$eig
	variance <- var_tbl.gen["constrained", "proportion"]
	p.val <- perm_anova.gen[1, 4]

	# extract the weighted average (sample) scores

	points <- capscale.gen$CCA$wa[, 1:2]
	points <- as.data.frame(points)
	colnames(points) <- c("x", "y")

	points <- cbind(points, design_temp[match(rownames(points), design_temp$SampleID_), ])

	# assign genotypes to the corresponding colors
	points$genotype    <- factor(points$genotype,    levels=genotype$names)

	# plot CPCo 1 and 2
	p <- ggplot(points, aes(x=x, y=y, colour=genotype, shape=rep)) +
		stat_ellipse(type="norm", linetype=2, aes(group=genotype), alpha=.3, level=.75) +
		geom_point(alpha=.7, size=2.5) +
		scale_colour_manual(values=genotype$colours, labels=geno.label) +
		labs(x=paste("CPCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
		     y=paste("CPCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
	     colour="Genotype", shape="Replicate") +
		ggtitle(paste(format(100 * variance, digits=3), " % of variance; p=",
			      format(p.val, digits=2),
			      sep="")) +
		theme_RTN +
		theme(legend.position="top",
		      legend.text.align=0)
	ggsave(p, file=paste(fig, "CPCoA_genotype.", ti, ".SP.pdf", sep=""), bg="transparent", width=4, height=4.5)



	k <- 2
	pcoa <- cmdscale(dist, k=k, eig=T)
	points <- pcoa$points
	eig <- pcoa$eig
	points <- as.data.frame(points)
	colnames(points) <- c("x", "y")

	points <- cbind(points, design_temp[match(rownames(points), design_temp$SampleID_), ])

	# assign genotypes to the corresponding colors
	points$genotype    <- factor(points$genotype,    levels=genotype$names)

	# plot CPCo 1 and 2
	p <- ggplot(points, aes(x=x, y=y, colour=genotype, shape=rep)) +
		stat_ellipse(type="norm", linetype=2, aes(group=genotype), alpha=.3, level=.75) +
		geom_point(alpha=.7, size=2.5) +
		scale_colour_manual(values=genotype$colours, labels=geno.label) +
		labs(x=paste("PC1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
		     y=paste("PC2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
	     colour="Genotype", shape="Replicate") +
		theme_RTN +
		theme(legend.position="top",
		      legend.text.align=0)
	ggsave(p, file=paste(fig, "PCoA_corZ.whole.", ti, ".SP.pdf", sep=""), bg="transparent", width=4, height=4.5)


}



# normal PCoA and PCoA constrained by timepoint for each genotype

g <- unique(design$genotype)
for(ge in g) {

	print(ge)

	idx <- design$genotype == ge
	design_temp <- design[idx, ]

	idx <- colnames(Z) %in% design_temp$SampleID_
	Z_temp <- Z[, idx]

	dist <- 1-cor(Z_temp)

	idx <- genotype$names == ge
	genotype_temp <- genotype[idx,]
	geno.label_temp <- geno.label[idx]


	sqrt_transform <- T
	d <- design_temp

	capscale.gen <- capscale(dist ~ timepoint + Condition(genotype + rep), data=d, add=F, sqrt.dist=sqrt_transform)

	# ANOVA-like permutation analysis
	perm_anova.gen <- anova.cca(capscale.gen)
	print(perm_anova.gen)
						    
	# generate variability tables and calculate confidence intervals for the variance
	var_tbl.gen <- variability_table(capscale.gen)
	eig <- capscale.gen$CCA$eig
	variance <- var_tbl.gen["constrained", "proportion"]
	p.val <- perm_anova.gen[1, 4]

	# extract the weighted average (sample) scores

	points <- capscale.gen$CCA$wa[, 1:2]
	points <- as.data.frame(points)
	colnames(points) <- c("x", "y")

	points <- cbind(points, design_temp[match(rownames(points), design_temp$SampleID_), ])

	# assign genotypes to the corresponding colors
	points$genotype    <- factor(points$genotype,    levels=genotype_temp$names)
	points$timepoint   <- factor(points$timepoint,   levels=timepoint$names)

	# plot CPCo 1 and 2
	p <- ggplot(points, aes(x=x, y=y, colour=genotype, shape=timepoint)) +
		stat_ellipse(type="norm", linetype=2, aes(group=timepoint), alpha=.3, level=.75) +
		geom_point(alpha=.7, size=2.5) +
		scale_colour_manual(values=genotype_temp$colours, labels=geno.label_temp, guide=F) +
		scale_shape_manual(values=timepoint$shapes) +
		labs(x=paste("CPCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
		     y=paste("CPCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
	     colour="Genotype", shape="Timepoint") +
		ggtitle(paste(format(100 * variance, digits=3), " % of variance; p=",
			      format(p.val, digits=2),
			      sep="")) +
		theme_RTN +
		theme(legend.position="top",
		      legend.text.align=0)
	ggsave(p, file=paste(fig, "CPCoA_timepoint.", ge, ".SP.pdf", sep=""), bg="transparent", width=4, height=4.5)




	k <- 2
	pcoa <- cmdscale(dist, k=k, eig=T)
	points <- pcoa$points
	eig <- pcoa$eig
	points <- as.data.frame(points)
	colnames(points) <- c("x", "y")

	points <- cbind(points, design_temp[match(rownames(points), design_temp$SampleID_), ])

	# assign genotypes to the corresponding colors
	points$genotype    <- factor(points$genotype,    levels=genotype_temp$names)
	points$timepoint   <- factor(points$timepoint,   levels=timepoint$names)

	# plot CPCo 1 and 2
	p <- ggplot(points, aes(x=x, y=y, colour=genotype, shape=timepoint)) +
		stat_ellipse(type="norm", linetype=2, aes(group=timepoint), alpha=.3, level=.75) +
		geom_point(alpha=.7, size=2.5) +
		scale_colour_manual(values=genotype_temp$colours, labels=geno.label_temp, guide=F) +
		scale_shape_manual(values=timepoint$shapes) +
		labs(x=paste("PC1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
		     y=paste("PC2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
	     colour="Genotype", shape="Timepoint") +
		theme_RTN +
		theme(legend.position="top",
		      legend.text.align=0)
	ggsave(p, file=paste(fig, "PCoA_corZ.whole.", ge, ".SP.pdf", sep=""), bg="transparent", width=4, height=4.5)


}








