## Scripts of Uemura and Nakano _et al._, A Golgi-released subpopulation of the trans-Golgi network mediates constitutive and pathogen-inducible protein secretion in Arabidopsis [Plant Physiol., 2019](http://www.plantphysiol.org/content/179/2/519).

originally by Ryohei Thomas Nakano

nakano@mpipz.mpg.de

These scripts are made available to facilitate reproducibility of our research. If you re-use any or part of this code, please reference with comments and cite our paper! Raw data and intermediate results necessary to run these scripts can be downloaded [here](http://www.mpipz.mpg.de/R_scripts).

---------------------------

### Scripts used for proteomic analysis
#### Data process

[1_proteinID.R](https://github.com/rtnakano1984/SYP4_secretome/blob/master/1_proteinID.R)  
R script for data loading and processing, merging all slpice variants. It exports a list of protein groups that are assigned with multiple loci. *The exported list of loci should be converted to a csv file by the GeneTable function at [ATTED-II](http://atted.jp/top_search.shtml#GeneTable)*

[2_GOcomp_proteinGroups.R](https://github.com/rtnakano1984/SYP4_secretome/blob/master/2_GOcomp_proteinGroups.R)  
R script for Data processing, dealing with missing loci in ATTED-II. *The exported list of protein groups should be manually curate, so that each protein group has single representative gene model.*

[3_proteinID_2.R](https://github.com/rtnakano1984/SYP4_secretome/blob/master/3_proteinID_2.R)  
R script for merging curated file with the original data file.

[4_data_process.R](https://github.com/rtnakano1984/SYP4_secretome/blob/master/4_data_process.R)  
R script for spliting the data file into files for each measures (LFQ, iBAQ, unique+razor peptides, and total peptides).  
*The exported list of loci is used to analyze their predicted signal paptides (SP) at [SignalP](http://www.cbs.dtu.dk/services/SignalP/), of which output is exported in the 'short' format*

[5_baseline_correction.R](https://github.com/rtnakano1984/SYP4_secretome/blob/master/5_baseline_correction.R)  
R script to calculate the bottom cutoff at 1% quantile and to replace any lower LFQ values by this value (Fig. S3B).

[6_data_process_SP.R](https://github.com/rtnakano1984/SYP4_secretome/blob/master/6_data_process_SP.R)  
R script for subseting the LFQ data table to those with predicted SPs.

[7_mat_Z.R](https://github.com/rtnakano1984/SYP4_secretome/blob/master/7_mat_Z.R)  
R script to calculate median-centered Z scores based on LFQ.  
　　
　　


#### Statistical analysis

[8_GLM.R](https://github.com/rtnakano1984/SYP4_secretome/blob/master/8_GLM.R)  
R script for fitting read counts to a generalized linear model.

[9_DEG_analysis.R](https://github.com/rtnakano1984/SYP4_secretome/blob/master/9_DEG_analysis.R)  
R script for extracting differentially abundant proteins and generate volcano plots (Fig. S3H).

[10_kmeanClust_AIC.SP.R](https://github.com/rtnakano1984/SYP4_secretome/blob/master/10_kmeanClust_AIC.SP.R)  
R script for deciding 'k' for the following k-means clustering based on AIC and BIC (Fig. S5A).

[11_kmeans.SP.R](https://github.com/rtnakano1984/SYP4_secretome/blob/master/11_kmeans.SP.R)  
R script for performing k-means clustering, followed by hierarchical clustering of the clusters (Fig. S5B).

[12_heatmaps.SP.R](https://github.com/rtnakano1984/SYP4_secretome/blob/master/12_heatmaps.SP.R)  
R script for plotting heamaps (Fig. 6A and S3I).

[13_kmeanClust_AIC.R](https://github.com/rtnakano1984/SYP4_secretome/blob/master/13_kmeanClust_AIC.R) / [14_kmeans.R](https://github.com/rtnakano1984/SYP4_secretome/blob/master/14_kmeans.R) / [15_heatmaps.R](https://github.com/rtnakano1984/SYP4_secretome/blob/master/15_heatmaps.R)  
Same as what Scripts #10-12 do but without subsetting to proteins with SPs (Fig. S3C-G).

[16_CPCoA.SP.R](https://github.com/rtnakano1984/SYP4_secretome/blob/master/16_CPCoA.SP.R)  
R script for Principal Coordinates Analysis (PCoA) and Canonical Analysis of Principal coodinates (CAP) (Fig. 6B,C and S4A-E).

[17_genotype_comparison.Z.SP.R](https://github.com/rtnakano1984/SYP4_secretome/blob/master/17_genotype_comparison.Z.SP.R)  
R script for comparing similarity between genotypes over the time course (Fig. S4F).

[18_topGO.SP.all.R](https://github.com/rtnakano1984/SYP4_secretome/blob/master/18_topGO.SP.all.R)  
R script for enrichment analysis of gene ontology (Fig. S5C,D).

[19_GeneSetEnrichmentAnalysis.R](https://github.com/rtnakano1984/SYP4_secretome/blob/master/19_GeneSetEnrichmentAnalysis.R) / [20_GeneSetEnrichmentAnalysis_plot.R](https://github.com/rtnakano1984/SYP4_secretome/blob/master/20_GeneSetEnrichmentAnalysis_plot.R)    
R script for Gene Set Enrichment Analysis (GSEA) and plotting its results as heatmaps (Fig. S6).
　　
　　


#### Other scripts

[plotting_parameters.R](https://github.com/rtnakano1984/129E_RNAseq/blob/master/plotting_parameters.R)  
R script that contains parameters required for the other scripts.

---------------------------

For any questions regarding these scripts, please contact

Ryohei Thomas Nakano

nakano@mpipz.mpg.de
