library(TimeSeriesAnalysis)
library(SummarizedExperiment)
library(ggplot2)

#Set up physical data in repository for example dataset
#This line of code can be removed if the example dataset is not used
write_example_data_to_dir(example_data='PBMC') #MURINE or PBMC

########### PARAM SET-UP ###########

#Give names to saved object and name of results repository
name_result_folder<-'TS_results_PBMC_example/'
obj_name<-'timeSeries_obj_PBMC.Rdata'

#Path to count data and sample data respectively
my_path_data<-'data/PBMC/raw_counts_TS'
my_path_sample_dta<-'data/PBMC/sample_file.csv'

#Set-up time series object parameters
diff_exp_type<-'DESeq2' #package used for DE – can also be 'limma'
p_val_filter_type<-'padj' #Either padj or pvalue, used to filter for significance
p_thresh<-0.05 #pvalue or padj value threshold for significance
l2fc_thresh<-1 #log(2)foldChange threshold for significance
name_control<-'LPS' #Name of experiment as seen in the sample file
name_experiment<-'IgM' #Name of control as seen in the sample file
graphic_vector<-c("#e31a1c","#1f78b4") #Pre-set colors for the groups

#Declare organism and load library
org_sem_sim='org.Hs.eg.db'
library('org.Hs.eg.db')

#Define specimen and ontology parameters
my_ont_gpro='GO:BP'
my_ont_sem_sim='BP'

#mmusculus
#hsapiens
#celegans
my_org_gpro='hsapiens' #Set the species for the gprofiler analysis

# The seed serves to create reproducible results with PART.
# A seed will ensure that the random components of PART clustering will be the same
# as long as the same seed is used. For more information on seeds, please consult
# this link: https://www.r-bloggers.com/2018/07/%F0%9F%8C%B1-setting-a-seed-in-r-when-using-parallel-simulation/
PART_seed=123456

PART_l2fc<-2 #log(2)foldChange threshold for PART clustering
PART_min_clust<-50 #Minimum cluster size for PART
PART_recursion<-100 #Number of recursions, default is 100

log_tp_traj<-FALSE #Defines if timepoints should be log transformed for illustration purposes

# Allows for all temporal combinations to be done instead
# of just sequential comparison. ex: do TP2vsTP1, TP3vsTP2, AND TP3vsTP1. In a normal instance
# only the first two comparison of the example would be run.
do_all_temporal_comparisons=FALSE

#The ancestors that will be queried, the ontology must be specified (BP,MF,or CC)
#Set to an empty vector c() if not required by the analysis
target_ancestors<-c('GO:0002253','GO:0019882','GO:0002404','GO:0002339','GO:0042386',
                    'GO:0035172','GO:0002252','GO:0006955','GO:0002520','GO:0090713',
                    'GO:0045321','GO:0001776','GO:0050900','GO:0031294','GO:0002262',
                    'GO:0002683','GO:0002684','GO:0002440','GO:0002682','GO:0002200',
                    'GO:0045058','GO:0002507')
ancestor_ontology<-'BP'

#Some extra set-up
name_save_obj<-paste0(name_result_folder,obj_name)#The object will be saved in result folder
#Create main directory for results
dir.create(name_result_folder)

my_group_names<-c(name_experiment,name_control)
names(graphic_vector)<-c(name_experiment,name_control)

########### END PARAMS ###########

#Some extra set-up
name_save_obj<-paste0(name_result_folder,obj_name)#The object will be saved in result folder
#Create main directory for results
dir.create(name_result_folder)

my_group_names<-c(name_experiment,name_control)
names(graphic_vector)<-c(name_experiment,name_control)



# OBJECT CREATION
if(obj_name %in% list.files(name_result_folder)==F){
  TS_object <- new('TimeSeries_Object',sample_data=prep_sample_data(my_path_sample_dta,my_group_names),
                   group_names=my_group_names,group_colors=graphic_vector,DE_method=diff_exp_type,
                   DE_p_filter=p_val_filter_type,DE_p_thresh=p_thresh,DE_l2fc_thresh=l2fc_thresh,
                   PART_l2fc_thresh=PART_l2fc,sem_sim_org=org_sem_sim,Gpro_org=my_org_gpro)
  # TS_object <- TS_load_example_data(TS_object,'PBMC')
  TS_object <- add_experiment_data(TS_object,sample_dta_path=my_path_sample_dta,count_dta_path=my_path_data)
  TS_object <- add_semantic_similarity_data(TS_object,my_ont_sem_sim)  TS_object <- add_semantic_similarity_data(TS_object,my_ont_sem_sim)
}else{
  load(name_save_obj)
}


#DE SECTION

#Perform normalization if the DESeq2 tool is being used and if normalized matrix doesn't exist
if (TS_object@DE_method=='DESeq2' & 'norm' %in% names(TS_object@count_matrix)!=T){
  TS_object <- normalize_timeSeries_with_deseq2(time_object=TS_object)
}

#Perform conditional differential gene expression analysis
TS_object<-conditional_DE_wrapper(TS_object)
#Perform temporal differential gene expression analysis
TS_object<-temporal_DE_wrapper(TS_object)

#save the timeseries object
save(TS_object,file=name_save_obj)


# PART SECTION
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
if (is.null(PART_seed)==FALSE){
  set.seed(as.character(PART_seed))
}

TS_object<-compute_PART(TS_object,part_recursion=PART_recursion,part_min_clust=PART_min_clust,
                        custom_seed=PART_seed,dist_param="euclidean", hclust_param="average")
#Save the TimeSeries object to the directory to prevent loss of data in the event
#of a downstream error
save(TS_object,file=name_save_obj)

#GPROFILER SECTION
# load('timeseries_obj_res.Rdata')
TS_object<-run_gprofiler_PART_clusters(TS_object) #Run the gprofiler analysis
#Save the results in a Rdata obbject
save(TS_object,file=name_save_obj)


