library(TimeSeriesAnalysis)
packages_for_loading<-c('DESeq2','limma','clusterGenomics','gprofiler2','tictoc',
                        'tibble','GOSemSim','dplyr','tibble')
lapply(packages_for_loading, require, character.only = TRUE)

#Set up physical data in repository for example dataset
#This line of code can be removed if the example dataset is not used
write_example_data_to_dir(example_data='PBMC') #MURINE or PBMC

########### PARAM SET-UP ###########

#Give names to saved object and name of results repository
name_result_folder<-'script_method/TS_results_example/'
obj_name<-'timeSeries_obj_example.Rdata'

# data
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

#Declare and install organism if necessary
org_sem_sim='org.Hs.eg.db'

my_org_gpro='hsapiens' #Set the species for the gprofiler analysis


PART_l2fc<-2 #log(2)foldChange threshold for PART clustering
PART_min_clust<-50 #Minimum cluster size for PART
PART_recursion<-100 #Number of recursions, default is 100, using 10 for example

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
  TS_object <- create_raw_count_matrix(TS_object,my_path_data)
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

#Use all samples, but implement a custom order. In this case it is reversed
samps_2<-TS_object@sample_data$sample[TS_object@sample_data$group==TS_object@group_names[2]]
samps_1<-TS_object@sample_data$sample[TS_object@sample_data$group==TS_object@group_names[1]]

#Create the matrix that will be used for PART clustering
TS_object<-prep_counts_for_PART(object=TS_object,target_genes=signi_genes,
                                scale=T,target_samples=c(samps_2,samps_1))

#Initiate PART clustering with a pre-defined seed – 123456
TS_object<-compute_PART(TS_object,part_recursion=PART_recursion,part_min_clust=PART_min_clust,
                        dist_param="euclidean", hclust_param="average",
                        custom_seed=123456)
#Save the TimeSeries object to the directory to prevent loss of data in the event
#of a downstream error
save(TS_object,file=name_save_obj)

#GPROFILER SECTION
TS_object<-run_gprofiler_PART_clusters(TS_object) #Run the gprofiler analysis
#Save the results in a Rdata obbject
save(TS_object,file=name_save_obj)


