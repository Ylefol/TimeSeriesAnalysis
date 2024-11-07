
#Load libraries and functions
library(TimeSeriesAnalysis)
library(ggplot2)
library(stringr)
library(DESeq2)
source('~/A_Projects/EpiGen/R_Work_Folder/TimeSeriesAnalysis/group_comparison_functions.R')


##### EXPLANATION
# This script retrieves differentially expressed genes from one or more instances
# of the TiSA pipeline. Note that this script only gives visual representation of multiple
# runs of the pipeline. Statistically speaking there are caveats, for example there will
# not be a correction for multiple group testing. The following figures serve as a means
# to compare the different runs of the pipeline with different settings or groups.
##### END

#Provide the path to several TiSA objects which will be used to retrieve differentially
#expressed genes.
TS_paths<-list('LPS-IgM'='rmarkdown_method/TS_results_PBMC_example/timeSeries_obj_PBMC_example.Rdata',
               'TGF-IgM'='rmarkdown_method/TS_results_PBMC_example_test_TGF/TiSA_object.Rdata')


#Provide group names which are also found in the objects given above.
#Give the group names in the order in which you would like them to appear
#in the heatmap. Also provide the color to be associated to each group
group_names<-c('IgM','LPS','TGF')
graphic_vector<-c("#1b9e77","#d95f02",'#7570b3')
names(graphic_vector)<-group_names

#Load up the sample data - provide the path
sample_data<-set_up_sample_dta('data/PBMC/sample_file.csv',group_names)

#This creates the deseq2 object and the normalized matrix. The matrix is normalized
#with all samples found in the sample_data (with it being filtered by group names)
path_to_counts<-'data/PBMC/raw_counts_TS'
returned_list<-create_matrix_and_dds(path_to_counts,sample_data)
dds<-returned_list[[1]]
norm_matrix<-returned_list[[2]]

#Retrieve the differential expressed genes
target_genes<-get_significant_genes(TS_object_paths=TS_paths,l2fc_thresh=1)

#Create a dummy TimeSeries object - Note that several parameters are given
#dummy or non-sensical values. This was done to illustrate that these will not be
#utilized in this script.
Dummy_obj <- new('TimeSeries_Object',
                 group_names=group_names,group_colors=graphic_vector,DE_method='DESeq2',
                 DE_p_filter='dummy',DE_p_thresh=999,DE_l2fc_thresh=999,
                 PART_l2fc_thresh=999,sem_sim_org='dummy',Gpro_org='hsapiens')
exp_data<-SummarizedExperiment(assays = norm_matrix,colData = sample_data)
names(assays(exp_data))='norm'
Dummy_obj@exp_data<-exp_data
Dummy_obj@DESeq2_obj<-dds

#Create the PCA plot
TS_pca<-plot_PCA_TS(Dummy_obj,DE_type='all')
TS_pca

#Prepare the subsetted counts with the differentially expressed genes. Load up the
#results in the Dummy object
Dummy_obj<-prep_counts_for_PART(object=Dummy_obj,target_genes=target_genes,
                                scale=TRUE,target_samples=sample_data$sample)

#Compute PART clustering and store in the Dummy object
Dummy_obj<-compute_PART(object=Dummy_obj,part_recursion=100,part_min_clust=50,
                          dist_param="euclidean", hclust_param="average",
                          custom_seed=123456, custom_matrix=NULL,return_as_object=TRUE)

##### NOTE
# All plots are made in such a way that the ggplot2 object is returned. This can
# be altered in some cases by specifying a path to which the results will be saved,
# however the easiest version would be to use ggsave() to save them in the desired
# height, width, dpi (resolution), and location.
##### END

#Plot the PART heatmap
#Set name to NULL to have the plot returned
PART_heat<-PART_heat_map(object=Dummy_obj,heat_name=NULL)
PART_heat


#Custom function to give more customization options - plots cluster trajectories
traj_plot<-create_cluster_traj_plot(object=Dummy_obj,
                                    target_clusters=c('C1','C2','C3','C4','C5','C6'),
                                    num_clusters_per_row=1,log_TP=FALSE)

traj_plot


#Plotting individual genes
gene_traj<-calculate_gene_traj_data(Dummy_obj,'AICDA',log_timepoint=FALSE)
aicda_plot<-plot_single_gene_traj(gene_traj,graphic_vector)
aicda_plot

#Run gprofiler
#Ensure that you have declared a valid organism in the dummy object (Gpro_org)
Dummy_obj<-run_gprofiler_PART_clusters(Dummy_obj)

#Set file_loc to NULL to return the plot instead of saving tp file
dotplot=GO_dotplot_wrapper(Dummy_obj,file_loc = NULL,target_ontology = 'REAC',top_n = 10,return_plot = TRUE)
dotplot


##### ANCESTOR QUERY
# The ancestor comparative query is a multi-faceted dotplot which shows the children
# of the queried ancestors that are found as significant in one or more of the
# provided TimeSeriesObjects. It also provides a directionality to each GO. Where an upward
# arrow shows that the cluster associated to the GO has a higher expression in the 'experiment'
# as opposed to the 'control'. The opposite would result in a downward arrow.

# If you would want to ignore the single timpoint parameter you can input 'mean' instead
# and thus the directionality of the arrow will be based on the 'mean' of the clusters
# values across timepoints.
##### END

ancestors_interest<-c('GO:0002253','GO:0019882','GO:0002404','GO:0002339','GO:0042386',
                    'GO:0035172','GO:0002252','GO:0006955','GO:0002520','GO:0090713',
                    'GO:0045321','GO:0001776','GO:0050900','GO:0031294','GO:0002262',
                    'GO:0002683','GO:0002684','GO:0002440','GO:0002682','GO:0002200',
                    'GO:0045058','GO:0002507')

# single timepoint
plot_title<-'Comparative dotplot for timepoint 1'
comp_ancestor_dta<-get_comparative_ancestor_data(TS_object_paths=TS_paths,target_ancestors=ancestors_interest,
                                                 target_timepoint=1,gpro_ont='GO:BP')
anc_plots<-plot_comp_ancestor_plots(plot_data=comp_ancestor_dta,title=plot_title)
anc_plots

# Mean of all timepoints
plot_title<-'Comparative dotplot for mean of timepoints'
comp_ancestor_dta<-get_comparative_ancestor_data(TS_object_paths=TS_paths,target_ancestors=ancestors_interest,
                                                 target_timepoint='mean',gpro_ont='GO:BP')
anc_plots<-plot_comp_ancestor_plots(plot_data=comp_ancestor_dta,title=plot_title)
anc_plots
