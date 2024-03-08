

#Note that if you utilize the results of this script, please cite the creators of clusterProfiler
#Citing information can be found on their github: https://github.com/YuLab-SMU/clusterProfiler


library(TimeSeriesAnalysis)
library(clusterProfiler)
library(ggplot2)
library(enrichplot)

#Load organism of choice
library('org.Hs.eg.db')

#Set target folder
target_folder<-'rmarkdown_method/TS_results_PBMC_example/'
load(paste0(target_folder,'timeSeries_obj_PBMC_example.Rdata'))

#Check possible names in conditional and temporal
names(TS_object@DE_results$conditional)
names(TS_object@DE_results$temporal)

#Select target DE results
target_DE='conditional'
target_name='IgM_vs_LPS_TP_3'

#Extract gene list
extracted_de <- TS_object@DE_results[[target_DE]][[target_name]]$DE_raw_data[,c('gene_id','log2FoldChange')]
gene_list<-extracted_de$log2FoldChange
names(gene_list)=extracted_de$gene_id
# omit any NA values
gene_list<-na.omit(gene_list)
# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)


#Check available keytypes
keytypes(org.Hs.eg.db)


#Run the GSEA. Important parameters are keyType, which is the keytype of
#the gene_list. In this case, SYMBOL. OrgDb must also be specified. Other
#parameters can be left as default or modified - more information is available
#on the cluster profiler page
gse <- gseGO(geneList=gene_list,
             ont ="ALL",
             keyType = "SYMBOL",
             nPerm = 10000,
             minGSSize = 3,
             maxGSSize = 800,
             pvalueCutoff = 0.05,
             verbose = TRUE,
             OrgDb = org.Hs.eg.db,
             pAdjustMethod = "none")
#Add the semantic similarity
gse<-pairwise_termsim(gse,semData=org.Hs.eg.db)

require(DOSE)
results_folder<-paste0(target_folder,'GSEA_results/',target_name,'/')
dir.create(results_folder,recursive = T)


#Below is a set of standard plots from clusterprofiler. These can be adjusted in many ways.
#First, the width and height should be adjusted on a per plot basis, as the required
#width and height will greatly vary depending on the results illustrated.

#The showCategory parameter is to show the top n results based on adjusted pvalue.

#You can also submit a vector of names to the showCategory parameter.
#This is illustrated in the last bit of code on this page

dotplot_name<-'default_dotplot.png'
dotplot<-dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
ggsave(plot = dotplot,filename = paste0(results_folder,dotplot_name),width = 9,height=10)

emap_name<-'default_emap.png'
emap<-emapplot(gse, showCategory = 10)
ggsave(plot = emap,filename = paste0(results_folder,emap_name),width = 7,height=7)

cnet_name<-'default_cnet.png'
cnet<-cnetplot(gse, categorySize="pvalue", color.params=list(foldChange=gene_list), showCategory = 3)
ggsave(plot = cnet,filename = paste0(results_folder,cnet_name),width = 10,height=10)

ridge_name<-'default_ridge.png'
ridge<-ridgeplot(gse) + labs(x = "enrichment distribution")
ggsave(plot = ridge,filename = paste0(results_folder,ridge_name),width = 9,height=15)



#Look at first 25 results
gse@result$Description[1:25]
#Creating a vector with three descriptions to illustrate use
desc_vector<-c(gse@result$Description[20:22])
#Create a cnet using desc_vector
cnet<-cnetplot(gse, categorySize="pvalue", color.params=list(foldChange=gene_list), showCategory = desc_vector)
cnet
