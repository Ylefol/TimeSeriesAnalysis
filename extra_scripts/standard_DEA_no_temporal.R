
setwd('~/A_Projects/EpiGen/R_Work_Folder/TimeSeriesAnalysis/')

library(TimeSeriesAnalysis)
#Give names to saved object and name of results repository
name_result_folder<-'TS_results_single_TP_test/'
obj_name<-'timeseries_obj_res.Rdata'

#Path to count data and sample data respectively
#Not needed if using example data
my_path_data<-'~/A_Projects/EpiGen/R_Work_Folder/TimeSeriesAnalysis/data/PBMC/raw_counts_TS/'
my_path_sample_dta<-'~/A_Projects/EpiGen/R_Work_Folder/TimeSeriesAnalysis/data/PBMC/sample_file_single_TP.csv'

org_sem_sim='org.Hs.eg.db'
library('org.Hs.eg.db')
name_control<-'LPS' #Name of control as seen in the sample file
name_experiment<-'IgM' #Name of experiment as seen in the sample file
graphic_vector<-c("#e31a1c","#1f78b4") #Pre-set colors for the groups
my_group_names<-c(name_experiment,name_control)
names(graphic_vector)<-c(name_experiment,name_control)

TS_object <- new('TimeSeries_Object',
                 group_names=my_group_names,group_colors=graphic_vector,DE_method='DESeq2',
                 DE_p_filter='padj',DE_p_thresh=0.05,DE_l2fc_thresh=1,
                 PART_l2fc_thresh=4,sem_sim_org=org_sem_sim,Gpro_org='hsapiens')

TS_object <- add_experiment_data(TS_object,sample_dta_path=my_path_sample_dta,count_dta_path=my_path_data)
TS_object <- add_semantic_similarity_data(TS_object,'BP')

TS_object <- normalize_timeSeries_with_deseq2(time_object=TS_object)
TS_object<-conditional_DE_wrapper(TS_object)
signi_genes<-select_genes_with_l2fc(TS_object)

sample_data<-exp_sample_data(TS_object)
#Use all samples, but implement a custom order. In this case it is reversed
samps_2<-sample_data$sample[sample_data$group==TS_object@group_names[2]]
samps_1<-sample_data$sample[sample_data$group==TS_object@group_names[1]]

#Create the matrix that will be used for PART clustering
TS_object<-prep_counts_for_PART(object=TS_object,target_genes=signi_genes,
                                scale=TRUE,target_samples=c(samps_2,samps_1))

#Sets a seed for reproducibility
PART_seed='123456'
if (is.null(PART_seed)==FALSE){
  set.seed(as.character(PART_seed))
}

TS_object<-compute_PART(TS_object,part_recursion=100,part_min_clust=50,
                        custom_seed=PART_seed,dist_param="euclidean", hclust_param="average")
TS_object<-run_gprofiler_PART_clusters(TS_object) #Run the gprofiler analysis
TS_pca<-plot_PCA_TS(TS_object,DE_type='all')

dir.create(name_result_folder)
plot_wrapper_DE_results(object=TS_object,DE_type='conditional',genes_of_interest=c(),results_folder=name_result_folder)
dir.create(paste0(name_result_folder,'/PART_results'))
PART_heat_map(TS_object,paste0(name_result_folder,'PART_results/PART_heat')) #Create a summary heatmap

#Create standard gprofiler results
gpro_res<-gprofiler_cluster_analysis(TS_object,'GO:BP',save_path = name_result_folder)
GO_clusters<-gpro_res[['GO_df']]
sem_dta<-slot(TS_object,'sem_list')
found_clusters<-find_clusters_from_termdist(GO_clusters,sem_dta)

#Plot and save MDS and clustered MDS
MDS_plots<-wrapper_MDS_and_MDS_clusters(GO_clusters,sem_dta,'BP',target_dir=paste0(name_result_folder,'gprofiler_results/'),return_plot=TRUE)


target_top=10
select_ontology<-'REAC'
GO_dotplot_wrapper(TS_object,file_loc=name_result_folder,target_ontology=select_ontology,top_n=target_top,custom_width=10)

target_ancestors<-c('GO:0002253','GO:0019882','GO:0002404','GO:0002339','GO:0042386',
                    'GO:0035172','GO:0002252','GO:0006955','GO:0002520','GO:0090713',
                    'GO:0045321','GO:0001776','GO:0050900','GO:0031294','GO:0002262',
                    'GO:0002683','GO:0002684','GO:0002440','GO:0002682','GO:0002200',
                    'GO:0045058','GO:0002507')
GOs_ancestors_clust<-find_relation_to_ancestors(target_ancestors,GO_clusters,ontology = 'BP')
ancestor_plots<-wrapper_ancestor_curation_plots(GOs_ancestors_clust,sem_dta,return_plot=TRUE,target_dir=name_result_folder)
