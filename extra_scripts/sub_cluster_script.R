library(TimeSeriesAnalysis)

library("org.Hs.eg.db")
my_org_sem_sim='org.Hs.eg.db'

log_tp_traj=F

targeted_cluster='C8'
ENSG_transcript_adjust=T

target_dir='rmarkdown_method/FOUNDIN-PD/TS_results_SNCA_Healthy/'

dir.create(paste0(target_dir,'/sub_cluster_analysis'))
sub_analysis_name=paste0(target_dir,'/sub_cluster_analysis/',targeted_cluster,'/')
dir.create(sub_analysis_name)

load(file = paste0(target_dir,'/timeseries_obj_res.Rdata'))

#Create the sub matrix
genes_of_cluster<-row.names(TS_object@PART_results$cluster_map)[TS_object@PART_results$cluster_map$cluster==targeted_cluster]
sub_matrix<-TS_object@PART_results$part_matrix[genes_of_cluster,]

#Re-do part for the sub matrix and replace necessary elements of the object
sub_TS_object<-compute_PART(TS_object,part_recursion=100,part_min_clust = 10,
                        dist_param="euclidean", hclust_param="average",
                        custom_seed=123456,custom_matrix=sub_matrix)

sub_TS_object@PART_results$part_matrix<-sub_matrix
sub_TS_object@Gprofiler_results<-list()
sub_TS_object<-run_gprofiler_PART_clusters(sub_TS_object,ENSG_transcript_adjust = ENSG_transcript_adjust)

save(sub_TS_object,file=paste0(sub_analysis_name,'/timeseries_obj_res.Rdata'))

load(file=paste0(sub_analysis_name,'/timeseries_obj_res.Rdata'))
#Create PART analysis related plots in a new sub-directory based on the above
#parameters
PART_heat_map(sub_TS_object,paste0(sub_analysis_name,'PART_heat'))

my_ts_data<-calculate_cluster_traj_data(sub_TS_object,scale_feat=T)
my_ts_data$group<-factor(my_ts_data$group,levels=rev(TS_object@group_names))
my_mean_ts_data<-calculate_mean_cluster_traj(my_ts_data)
my_mean_ts_data$group<-factor(my_mean_ts_data$group,levels=rev(TS_object@group_names))

wrapper_cluster_trajectory(sub_TS_object,my_ts_data,my_mean_ts_data,log_TP=log_tp_traj,plot_name=paste0(sub_analysis_name,'Ctraj'))

#Params
my_ont_gpro='GO:BP'
clusts_int<-unique(sub_TS_object@PART_results$cluster_map$cluster)
GO_clusters<-gprofiler_cluster_analysis(sub_TS_object,my_ont_gpro,save_path = sub_analysis_name)
GO_clusters<-GO_clusters$GO_df
#Filter for modules of interest
GO_clusters<-GO_clusters[GO_clusters$group_name %in% clusts_int,]

#Create semantic data
sem_data <- slot(TS_object,'sem_list')

#Plot and save MDS and clustered MDS
wrapper_MDS_and_MDS_clusters(GO_clusters,sem_data,'BP',target_dir=paste0(sub_analysis_name,'gprofiler_results/'))

#Create a dotplot summarizing the top n findings for each cluster
target_top=10
select_ontologies<-c('GO:MF','GO:BP','GO:CC','REAC','CORUM','HPA','KEGG','WP','TF','MIRNA','HP','WP')
for(ont in select_ontologies){
  GO_dotplot_wrapper(sub_TS_object,file_loc=sub_analysis_name,target_ontology=ont,top_n=target_top,custom_width=10)
}


# #Dotplot for terms relating to specific ancestors
# #All GOs relating to immune system process
# target_ancestors<-c('GO:0002253','GO:0019882','GO:0002404','GO:0002339','GO:0042386',
#                     'GO:0035172','GO:0002252','GO:0006955','GO:0002520','GO:0090713',
#                     'GO:0045321','GO:0001776','GO:0050900','GO:0031294','GO:0002262',
#                     'GO:0002683','GO:0002684','GO:0002440','GO:0002682','GO:0002200',
#                     'GO:0045058','GO:0002507')
# #sensory perception of small and G protein coupled receptor signalling pathway
# # target_ancestors<-c('GO:0007608','GO:0007186')
#
# GOs_ancestors_clust<-find_relation_to_ancestors(target_ancestors,GO_clusters,ontology = 'BP')
# wrapper_ancestor_curation_plots(GOs_ancestors_clust,sem_data,target_dir=sub_analysis_name)
