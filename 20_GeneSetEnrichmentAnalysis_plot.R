
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
library(org.At.tair.db, quietly=T, warn.conflicts=F)
library(GOSemSim,       quietly=T, warn.conflicts=F)
library(stringr,        quietly=T, warn.conflicts=F)
library(dplyr,          quietly=T, warn.conflicts=F)
library(pipeR,          quietly=T, warn.conflicts=F)
library(reshape2,       quietly=T, warn.conflicts=F)
library(ggplot2,        quietly=T, warn.conflicts=F)

# directories
original_data  <- "/your/data/directory/original_data/"
processed_data <- "/your/data/directory/processed_data/"
fig            <- "/your/data/directory/fig/"
stat           <- "/your/data/directory/statistics/"
scripts        <- "/your/data/directory/scripts/"
source(paste(scripts, "plotting_parameters.R", sep=""))


# data import
design       <-           read.table(paste(processed_data, "design.txt", sep=""),   sep="\t", header=T, check.names=F, stringsAsFactors=F)
gse          <-           read.table(paste(stat, "GSEO_all.txt", sep=""),           sep="\t", header=T, check.names=F, stringsAsFactors=F)
logFC_P      <- as.matrix(read.table(paste(processed_data, "logFC.P.txt", sep=""),  sep="\t", header=T, check.names=F, stringsAsFactors=F, row.names=1))
SP_genes     <-           read.table(paste(processed_data, "SP_genes.txt", sep=""), sep="\t", header=T, check.names=F, stringsAsFactors=F)[,1]

# merge
idx <- match(gse$contrasts, design$contrast)
gse <- data.frame(gse, design[idx, ], stringsAsFactors=F)

terms <- unique(gse[, c("ID", "ONTOLOGY", "Description")])
terms$Description <- as.character(terms$Description)


# sort by semantic similarity of GO terms

sorted_ID_list <- lapply(c("BP", "CC", "MF"), function(x) {
	db <- godata("org.At.tair.db", ont=x)
	term <- terms$ID[terms$ONTOLOGY == x]
	sim <- mgoSim(term, term, semData=db, measure="Wang", combine=NULL)
	d <- as.dist(1-sim)
	hclust <- hclust(d, "ave")
	sorted <- rownames(sim)[hclust$order]
	return(sorted)
})

sorted_simGO <- do.call(c, sorted_ID_list)

idx <- match(sorted_simGO, terms$ID)
sorted_simGO_des <- terms$Description[idx]



# sort

# sort
gse$SampleID    <- factor(gse$SampleID,                  levels=design$SampleID)
gse$Description <- factor(as.character(gse$Description), levels=sorted_simGO_des)
gse$ID          <- factor(as.character(gse$ID),          levels=sorted_simGO)
gse$timepoint   <- factor(gse$timepoint,                 levels=timepoint$names)
gse$genotype    <- factor(gse$genotype,                  levels=genotype$names)

idx <- gse$p.adjust < alpha
gse$sig[!idx] <- 0
gse$sig[idx]  <- 1
gse$sig <- factor(gse$sig, levels=c(0,1))



# heatmap of NES
p <- ggplot(gse, aes(x=SampleID, y=Description, fill=NES)) +
	geom_raster() +
	scale_fill_gradient2(high=c_dark_green, mid=c_white, low=c_cudo_magenta, midpoint=0, na.value="gray95") +
	labs(x="", y="", fill="Normalized Enrichment Score") +
	facet_grid(ONTOLOGY ~ genotype, space="free", scales="free", drop=T, labeller=label_parsed) +
	theme_RTN +
	theme(legend.position="top",
		strip.text.x=element_blank(),
		axis.text.x=element_text(size=6, angle=30, hjust=1, vjust=1),
		axis.text.y=element_text(size=2.5),
		axis.line=element_blank(),
		panel.background=element_rect(colour=c_black, size=1, fill=NA),
		axis.ticks=element_blank()) +
	guides(fill=guide_colourbar(barheight=1, barwidth=5))
ggsave(p, file=paste(fig, "GSE_nomalizedES.sorted_simGO.png", sep=""), bg="transparent", width=6, height=20)

p <- p + scale_fill_gradient2(low=c_blue, mid=c_white, high=c_red, midpoint=0, na.value="grey95")
ggsave(p, file=paste(fig, "GSE_nomalizedES.sorted_simGO.blue_red.png", sep=""), bg="transparent", width=6, height=20)


# heatmap of p
p <- ggplot(gse, aes(x=SampleID, y=Description, fill=factor(sig))) +
	geom_raster() +
	scale_fill_manual(values=c("0"=c_white, "1"=c_red)) +
	labs(x="", y="", fill=paste("Significance (P < ", alpha, sep="")) +
	facet_grid(ONTOLOGY ~ genotype, space="free", scales="free", drop=T, labeller=label_parsed) +
	theme_RTN +
	theme(legend.position="top",
		strip.text.x=element_blank(),
		axis.text.x=element_text(size=6, angle=30, hjust=1, vjust=1),
		axis.text.y=element_text(size=2.5),
		axis.line=element_blank(),
		panel.background=element_rect(colour=c_black, size=1, fill=NA),
		axis.ticks=element_blank()) +
	guides(fill=guide_colourbar(barheight=1, barwidth=5))
ggsave(p, file=paste(fig, "GSE_padjust.sorted_simGO.png", sep=""), bg="transparent", width=6, height=20)



# plot heatmap of FC
idx <- rownames(logFC_P) %in% SP_genes
logFC_P <- logFC_P[idx, ]

ID_table <- AnnotationDbi::select(org.At.tair.db, keys=sorted_simGO, keytype="GOALL", columns="TAIR")

idx <- is.na(ID_table$TAIR)
ID_table <- ID_table[!idx,]

logFC_list <- lapply(unique(ID_table$GO), function(x) {
		message(x)
		idx <- ID_table$GO == x
		IDs <- ID_table$TAIR[idx]
        
        idx <- rownames(logFC_P) %in% IDs
        logFC <- logFC_P[idx, str_detect(colnames(logFC_P), "logFC")]
        if(length(IDs) == 1) {
                melt <- data.frame(Var1=IDs, Var2=names(logFC), value=logFC, Term=x)
        } else {
                melt <- melt(logFC)
                melt$Term <- x
        }

        melt$Var2 <- str_replace(melt$Var2, "_logFC", "")

        return(melt)
        })

logFC <- rbind_all(logFC_list)
logFC_ag <- logFC %>>% group_by(Term, Var2) %>>% summarise(logFC=mean(value))


# merge
idx <- match(logFC_ag$Var2, design$contrast)
logFC_ag <- data.frame(logFC_ag, design[idx, ], stringsAsFactors=F)

idx <- match(logFC_ag$Term, terms$ID)
logFC_ag$Description <- terms$Description[idx]
logFC_ag$ONTOLOGY    <- terms$ONTOLOGY[idx]


# sort
logFC_ag$Description <- factor(logFC_ag$Description, levels=sorted_simGO_des)
logFC_ag$Term        <- factor(logFC_ag$Term,        levels=sorted_simGO)
logFC_ag$SampleID    <- factor(logFC_ag$SampleID,  levels=design$SampleID)
logFC_ag$timepoint   <- factor(logFC_ag$timepoint, levels=timepoint$names)
logFC_ag$genotype    <- factor(logFC_ag$genotype,  levels=genotype$names)
levels(logFC_ag$genotype) <- geno.label


# saturate logFC
max <- quantile(logFC_ag$logFC, .99)
idx <- (logFC_ag$logFC > max)
logFC_ag$logFC[idx] <- max

min <- quantile(logFC_ag$logFC, .01)
idx <- (logFC_ag$logFC < min)
logFC_ag$logFC[idx] <- min

p <- ggplot(logFC_ag, aes(x=SampleID, y=Description, fill=logFC)) +
        geom_raster() +
        scale_fill_gradient2(high=c_dark_green, mid="white", low=c_cudo_magenta, midpoint=0, na.value="gray95") +
        labs(x="", y="", fill="logFC") +
        facet_grid(ONTOLOGY ~ genotype, space="free", scales="free", drop=T, labeller=label_parsed) +
        theme_RTN +
	theme(legend.position="top",
		strip.text.x=element_blank(),
		axis.text.x=element_text(size=6, angle=30, hjust=1, vjust=1),
		axis.text.y=element_text(size=2.5),
		axis.line=element_blank(),
		panel.background=element_rect(colour=c_black, size=1, fill=NA),
		axis.ticks=element_blank()) +
	guides(fill=guide_colourbar(barheight=1, barwidth=5))
ggsave(p, file=paste(fig, "GSE_logFC.sorted_simGO.png", sep=""), bg="transparent", width=6, height=20)

p <- p + scale_fill_gradient2(low=c_blue, mid=c_white, high=c_red, midpoint=0, na.value="grey95")
ggsave(p, file=paste(fig, "GSE_logFC.sorted_simGO.blue_red.png", sep=""), bg="transparent", width=6, height=20)

