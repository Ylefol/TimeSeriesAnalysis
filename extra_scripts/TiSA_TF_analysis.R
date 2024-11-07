
library(TimeSeriesAnalysis)
setwd('~/A_Projects/EpiGen/R_Work_Folder/TimeSeriesAnalysis/')
load('results_folder/TS_results_PBMC_example_test_script_method/TiSA_object.Rdata')

source('extra_scripts/TiSA_TF_functions.R')

#use 'human' or 'mouse'
organism='human'
save_location='results_folder/TS_results_PBMC_example_test_script_method/TF_results/'
dir.create(save_location)

#Gather TF database
my_net <- get_collectri(organism=organism, split_complexes=FALSE)

#Basic run of decouplR using either counts or DEGs to find significant transcription factors - top 25
my_found_TFs<-wrapper_TF_extraction(time_object = TS_object,TF_network = my_net,TF_method = 'DEG',save_location=save_location,top_TF=25)
my_found_TFs_count<-wrapper_TF_extraction(time_object = TS_object,TF_network = my_net,TF_method = 'counts',save_location=save_location,top_TF=25)

#Load JASPAR 2024 dataset
filename <- "data/JASPAR2024.sqlite"
sqlite.driver <- dbDriver("SQLite")
db <- dbConnect(sqlite.driver,
                dbname = filename)

#Set the identifier for the species in the JASPAR database
if(organism=='human'){
  target_species_id<-9606
}else if(organism=='mouse'){
  target_species_id<-10090
}


#Add nucleotide content - we use AT
AT_table<-calculate_nucleotide_content(JASPAR_dta=db,species_ID=target_species_id,
                                       target_nucleotide=c('A','T'))

#Run the decouplR analysis to add pvalues to all TFs of the table
returned_list<-add_pvalue_from_decouplr(time_object=TS_object,TF_network=my_net,nuc_table=AT_table)
AT_table<-returned_list[['nucleotide_table']]
TF_full_df<-returned_list[['TF_full_df']]

#Select significant TFs with the highest nucleotide content.
target_TFs_deg<-unique(AT_table$NAME)[0:100]
TF_full_df<-TF_full_df[TF_full_df$source %in% target_TFs_deg,]

#Add the nucleotide content/percentage to the names of the TF - gives better visualization
TF_with_AT<-add_nucleotide_percentage_to_TF(AT_table,TF_full_df)

#Plot the results
DE_list<-create_DE_list_from_time_object(TS_object)
create_deg_TF_heatmap(DE_list,TF_with_AT,custom_group_names,custom_colors,
                      save_name=paste0(save_location,'DE_heatmap_with_content.pdf'))


#Plotting volcano plot for highest content TF
target_TF<-unique(AT_table$NAME)[1]
TF_volcano<-plot_individual_TF(TS_object,target_TF=target_TF,net=my_net,DE_type='conditional',target_exp='IgM_vs_LPS_TP_1',save_location=save_location)


#Need to handle the sequence plotting
target_TF<-unique(AT_table$ID)[1]
plot_seqlogo(JASPAR_dta=db,target_TF_ID=target_TF,save_location=save_location)


