
library(TimeSeriesAnalysis)
library(WGCNA)

source('extra_scripts/WGCNA/WGCNA_custom_functions.R')

#Load the time object containing the clusters of interest
load('rmarkdown_method/TS_results_CoV_Mech/timeSeries_obj_CoV_example.Rdata')

#Load the metadata
datTraits<-read.csv('data/TS_covid/cov_mech_metadata.csv',row.names=1)

#Get gene matrix (normalized)
datExpr0<-t(assays(TS_object@exp_data)$norm)
sampleTree = hclust(dist(datExpr0), method = "average")
# Determine cluster under the line, re-cluster if necessary
returned_list <- cluster_under_line(sampleTree,datExpr0,datTraits)
#Assign the returned values to variables - only need expression matrix
datExpr <- returned_list[[3]]

#Plot the power plot
plot_power(datExpr)

#Manually set the power based on the above plot
#NOTE: The power is set by taking the first number that 'plateaus' above the red line
target_power<-8

#ensure the matrix has the right format
datExpr_comp <- matrix(as.numeric(datExpr),ncol = ncol(datExpr))#Ensure that matrix contains numerics

#Perform the WGCNA block clustering
#The below function will likely take 2 to 3 hours to run depending on the dataset
net = blockwiseModules(datExpr_comp, power = target_power,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = F,
                       maxBlockSize = 35000,
                       verbose = 3,nThreads=2)

save(net,file = 'extra_scripts/WGCNA/blockWiseModules.Rdata')

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;

#Perform correlation using pearson correlation
nSamples = nrow(datExpr)

MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)


#Plot the heatmap
png("labelled_heat_map_complete.png",25,15,units='cm',res=600)
plot_labeled_heatmap(moduleTraitCor,moduleTraitPvalue,datTraits,MEs)
dev.off()



