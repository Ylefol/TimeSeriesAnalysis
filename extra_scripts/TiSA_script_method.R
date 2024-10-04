
#Load libraries
library(TimeSeriesAnalysis)
library(SummarizedExperiment)
library(ggplot2)

#This writes example data to directory, enable it to run the example below.
#Only needs to be run once
#write_example_data_to_dir('PBMC')

#Set the name where the results will be stored and create dir
name_result_folder<-'TS_results_PBMC_example_test_script_method/'
dir.create(name_result_folder)

#Set graphics vector for group coloring
graphic_vector<-c("#e31a1c","#1f78b4") #Pre-set colors for the groups
names(graphic_vector)<-c('IgM','LPS')

#Declare organism and load library, install if necessary
library('org.Hs.eg.db')

# Object creation -------------------
TS_object <- new('TimeSeries_Object',
                 group_names=c('IgM','LPS'),group_colors=graphic_vector,DE_method='DESeq2',
                 DE_p_filter='padj',DE_p_thresh=0.05,DE_l2fc_thresh=1,
                 PART_l2fc_thresh=4,sem_sim_org='org.Hs.eg.db',Gpro_org='hsapiens')

TS_object <- add_experiment_data(TS_object,sample_dta_path='data/PBMC/sample_file.csv',count_dta_path='data/PBMC/raw_counts_TS')
TS_object <- add_semantic_similarity_data(TS_object,'BP')


# DESeq2 -------------------
TS_object <- normalize_timeSeries_with_deseq2(time_object=TS_object)
#Perform conditional differential gene expression analysis
TS_object<-conditional_DE_wrapper(TS_object)
#Perform temporal differential gene expression analysis
TS_object<-temporal_DE_wrapper(TS_object,do_all_combinations=FALSE)

# PART -------------------
#Extract genes for PART clustering based on defined log(2)foldChange threshold
signi_genes<-select_genes_with_l2fc(TS_object)

sample_data<-exp_sample_data(TS_object)
#Use all samples, but implement a custom order. In this case it is reversed
samps_2<-sample_data$sample[sample_data$group==TS_object@group_names[2]]
samps_1<-sample_data$sample[sample_data$group==TS_object@group_names[1]]

#Create the matrix that will be used for PART clustering
TS_object<-prep_counts_for_PART(object=TS_object,target_genes=signi_genes,
                                scale=TRUE,target_samples=c(samps_2,samps_1))

#Sets a seed for reproducibility
set.seed('123456')


TS_object<-compute_PART(TS_object,part_recursion=100,part_min_clust=50,
                        custom_seed=123456,dist_param="euclidean", hclust_param="average")

# Gprofiler -------------------
TS_object<-run_gprofiler_PART_clusters(TS_object) #Run the gprofiler analysis

# PCA -------------------
TS_pca<-plot_PCA_TS(TS_object,DE_type='all')
ggsave(paste0(name_result_folder,"PCA_plot.png"),dpi=300,width=21, height=19, units='cm',plot=TS_pca)

# DESeq2 results -------------------
#Set genes of interest (optional) - can be left as c()
genes_of_interest <- c('AICDA','APOBEC3H','APOBEC3F','APOBEC3D','APOBEC3C','APOBEC3G','APOBEC3B','APOBEC3A','SMUG1','UNG','EGFR')
#Run wrappers twice once for conditional and another for temporal
plot_wrapper_DE_results(object=TS_object,DE_type='conditional',genes_of_interest=genes_of_interest,results_folder=name_result_folder)
plot_wrapper_DE_results(object=TS_object,DE_type='temporal',genes_of_interest=genes_of_interest,results_folder=name_result_folder)

# PART results -------------------
dir.create(paste0(name_result_folder,'PART_results')) #create the directory to store results
PART_heat_map(TS_object,paste0(name_result_folder,'PART_results/PART_heat')) #Create a summary heatmap
ts_data<-calculate_cluster_traj_data(TS_object,scale_feat=TRUE) #Calculate scaled gene values for genes of clusters
mean_ts_data<-calculate_mean_cluster_traj(ts_data) #Calculate the mean scaled values for each cluster

#Function which determines the number of SVG files to plot all cluster trajectories
wrapper_cluster_trajectory(TS_object,ts_data,mean_ts_data,yaxis_name='scaled mean counts',log_TP=F,plot_name=paste0(name_result_folder,'/PART_results/Ctraj'))

# Gprofiler results -------------------
#Create standard gprofiler results
gpro_res<-gprofiler_cluster_analysis(TS_object,'GO:BP',save_path = name_result_folder)
GO_clusters<-gpro_res[['GO_df']]
sem_dta<-slot(TS_object,'sem_list')
found_clusters<-find_clusters_from_termdist(GO_clusters,sem_dta)

#Plot and save MDS and clustered MDS
MDS_plots<-wrapper_MDS_and_MDS_clusters(GO_clusters,sem_dta,'BP',target_dir=paste0(name_result_folder,'gprofiler_results/'),return_plot=TRUE)

#Create dotplots
GO_dotplot_wrapper(TS_object,file_loc=name_result_folder,target_ontology='GO:BP',top_n=10,custom_width=10)


# Ancestor queries results -------------------
GOs_ancestors_clust<-find_relation_to_ancestors(target_ancestors=c('GO:0002253','GO:0019882'),GOs_to_check=GO_clusters,ontology = 'BP')
ancestor_plots<-wrapper_ancestor_curation_plots(GOs_ancestors_clust,sem_dta,return_plot=TRUE,target_dir=name_result_folder)

save(TS_object,file=paste0(name_result_folder,'TiSA_object.Rdata'))

