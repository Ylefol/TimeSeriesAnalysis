
library(TimeSeriesAnalysis)
library(WGCNA)

source('extra_scripts/WGCNA/WGCNA_custom_functions.R')

#Load the time object containing the clusters of interest
load('rmarkdown_method/TS_results_CoV_Mech/timeSeries_obj_CoV_example.Rdata')

#Load the metadata
datTraits<-read.csv('data/TS_covid/cov_mech_metadata.csv',row.names=1)

#Some matrix setup
TS_part_genes<-row.names(TS_object@PART_results$part_data)
TS_corr_matrix<-t(assays(TS_object@exp_data)$norm[TS_part_genes,])
TS_corr_vector<-TS_object@PART_results$part_data$gene_cluster
nSamples = nrow(TS_corr_matrix)

#Get correlation based on first Principal component of each cluster - uses pearsons correlation
MEs0_custom = moduleEigengenes(TS_corr_matrix, TS_corr_vector)$eigengenes
MEs_custom = orderMEs(MEs0_custom)
moduleTraitCor_custom = cor(MEs_custom, datTraits, use = "p")
moduleTraitPvalue_custom = corPvalueStudent(moduleTraitCor_custom, nSamples)

#Plots and saves the labeled heatmap
png("labelled_heat_map_PART.png",25,15,units='cm',res=600)
plot_labeled_heatmap(moduleTraitCor_custom,moduleTraitPvalue_custom,datTraits,MEs_custom)
dev.off()


