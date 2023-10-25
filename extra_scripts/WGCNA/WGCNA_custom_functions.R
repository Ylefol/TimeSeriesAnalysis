plot_labeled_heatmap <- function(moduleTraitCor,moduleTraitPvalue,datTraits,MEs){
  # # Will display correlations and their p-values
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(6, 8.5, 3, 3));
  # Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(datTraits),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colorMatrix = NULL,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))

}




plot_power <-function(datExpr){
  #power settings
  powers = c(seq(1, 10, by = 1), seq(12, 30, by = 2))
  # Call the network topology analysis function
  sft = pickSoftThreshold(datExpr,power=powers,verbose=5)#, power=powers, verbose = 5)

  # sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")
  abline(h=0.80,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
}


cluster_under_line <- function(sample_tree,datExpr0,datTraits,cutoff_inputted=NULL,cluster_chosen=0){

  if (is.null(cutoff_inputted)==F){
    clust = cutreeStatic(sampleTree, cutHeight = cutoff_inputted, minSize = 10)
  }
  clust = cutreeStatic(sampleTree, minSize = 10)
  # table(clust)

  keepSamples = (clust==cluster_chosen)
  datExpr = datExpr0[keepSamples, ]
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)

  # Re-cluster samples
  sampleTree2 = hclust(dist(datExpr), method = "average")
  # Convert traits to a color representation: white means low, red means high, grey means missing entry
  traitColors = numbers2colors(datTraits, signed = FALSE)

  return_list<-list(sampleTree2,traitColors,datExpr)
  return(return_list)
}
