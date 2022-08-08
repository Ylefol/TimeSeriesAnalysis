library(TimeSeriesAnalysis)
packages_for_loading<-c('DESeq2','limma','grid','ggrepel','ggplot2','ComplexHeatmap',
                        'reshape2','plotly','GO.db','GOSemSim','data.table','gprofiler2',
                        'htmlwidgets','dynamicTreeCut','stringr')

lapply(packages_for_loading, require, character.only = TRUE)

#If timepoints should be log transformed for plotting
log_tp_traj=F

#Give names to saved object and name of results repository
name_result_folder<-'script_method/TS_results_example/'
obj_name<-'timeSeries_obj_example.Rdata'

load(file = paste0(name_result_folder,obj_name))

#Calcualte and plot the TimeSeries PCA
TS_pca<-plot_PCA_TS(TS_object,DE_type='all')
ggsave(paste0(name_result_folder,"PCA_plot.png"),dpi=300,width=21, height=19, units='cm',plot=TS_pca)

#Create a list of genes of interest
genes_of_interest <- c('AICDA','APOBEC3H','APOBEC3F','APOBEC3D','APOBEC3C','APOBEC3G','APOBEC3B','APOBEC3A','SMUG1','UNG','EGFR')

#Create results for conditional and temporal differential gene expression results
plot_wrapper_DE_results(object=TS_object,DE_type='conditional',genes_of_interest=genes_of_interest,results_folder=name_result_folder)
plot_wrapper_DE_results(object=TS_object,DE_type='temporal',genes_of_interest=genes_of_interest,results_folder=name_result_folder)

#Create the repository for genes_of_interest and perform the analysis
dir.create(paste0(name_result_folder,'genes_of_interest/'))
#May create warnings if there are too few timepoints, plots will be created as intended
create_tables_genes_of_interest_DE(TS_object,genes_of_interest,save_location=paste0(name_result_folder,'genes_of_interest/'),log_tp = log_tp_traj)

dir.create(paste0(name_result_folder,'PART_results')) #create the directory to store results
PART_heat_map(TS_object,paste0(name_result_folder,'PART_results/PART_heat')) #Create a summary heatmap
ts_data<-calculate_cluster_traj_data(TS_object,scale_feat=T) #Calculate scaled gene values for genes of clusters
mean_ts_data<-calculate_mean_cluster_traj(ts_data) #Calculate the mean scaled values for each cluster

#Function which determines the number of SVG files to plot all cluster trajectories
wrapper_cluster_trajectory(TS_object,ts_data,mean_ts_data,log_TP=log_tp_traj,plot_name=paste0(name_result_folder,'/PART_results/Ctraj'))



##########################
# G profiler results
##########################
#Define specimen and ontology parameters
library(TS_object@sem_sim_org,character.only = T)
my_ont_sem_sim='BP'
my_ont_gpro=paste0('GO:',my_ont_sem_sim)


#Create standard gprofiler results
gpro_res<-gprofiler_cluster_analysis(TS_object,my_ont_gpro,save_path = name_result_folder)
GO_clusters<-gpro_res[['GO_df']]

#Create semantic data
sem_data <- godata(TS_object@sem_sim_org, ont=my_ont_sem_sim, computeIC=TRUE)

#Plot and save MDS and clustered MDS
wrapper_MDS_and_MDS_clusters(GO_clusters,sem_data,my_ont_sem_sim,target_dir=paste0(name_result_folder,'gprofiler_results/'))



#Create a dotplot summarizing the top n findings for each cluster
target_top=10

select_ontology<-'GO:BP'
GO_dotplot_wrapper(gpro_file_location=name_result_folder,target_ontology=select_ontology,top_n=target_top)

select_ontology<-'REAC'
GO_dotplot_wrapper(gpro_file_location=name_result_folder,target_ontology=select_ontology,top_n=target_top)

select_ontology<-'KEGG'
GO_dotplot_wrapper(gpro_file_location=name_result_folder,target_ontology=select_ontology,top_n=target_top)

#Dotplot for terms relating to specific ancestors
#All GOs relating to immune system process
target_ancestors<-c('GO:0002253','GO:0019882','GO:0002404','GO:0002339','GO:0042386',
                    'GO:0035172','GO:0002252','GO:0006955','GO:0002520','GO:0090713',
                    'GO:0045321','GO:0001776','GO:0050900','GO:0031294','GO:0002262',
                    'GO:0002683','GO:0002684','GO:0002440','GO:0002682','GO:0002200',
                    'GO:0045058','GO:0002507')


ancestor_ontology<-'BP'

GOs_ancestors_clust<-find_relation_to_ancestors(target_ancestors,GO_clusters,ontology = ancestor_ontology)
ancestor_plots<-wrapper_ancestor_curation_plots(GOs_ancestors_clust,sem_data,target_dir=name_result_folder)

