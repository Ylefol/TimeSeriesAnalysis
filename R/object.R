# library(DESeq2)
# # library(limma)
# packages_for_loading<-c('DESeq2','limma')
# suppressPackageStartupMessages(lapply(packages_for_loading, require, character.only = TRUE))
#

# Class creation  -------------------


#' @title The TimeSeries class
#'
#' @description The TimeSeries class is the main class for the timeseries pipeline
#' It is used to store the input data along with the final processed data and
#' any intermediate results that may be used for plot creation
#'
#' @slot count_matrix A list used to store the raw and normalized count matrices
#' @slot sample_data A dataframe created by the user indicating the grouping, timepoints,
#' and replicates to which the samples belong
#' @slot group_names A vector for the group names used, should be in order for
#' which the differential expression comparisons will occur, where the first
#' should be the experiment and second should be the control
#' @slot group_colors A named vector indicating which colors should be associated to which group
#' @slot DE_method Either limma or DESeq2 to indicate which differential expression method is used
#' @slot DE_p_filter Character of paj or pvalue to use in the establishment of
#' significant differential expressed genes
#' @slot DE_p_thresh A numeric value used along with the pvalue filter to establish
#' significance. Genes will be considered significant if they are below this threshold
#' @slot DE_l2fc_thresh A numeric value for the log2FolChange threshold to establish
#' significance where the absolute value of the genes must be greater than the threshold
#' @slot Gpro_org The organism used in gprofiler friendly format
#' @slot sem_list The list containing the semantic similarity measures for each ontology
#' @slot DESeq2_obj The normalized DESeq2 object
#' @slot limma_object The EList limma object
#' @slot DE_results A list of results for the different differential expression experiments
#' performed
#' @slot PART_L2FC_thresh A integer indicating the log(2)foldchange threshold for genes
#' to be PART clustered
#' @slot PART_results A list of results for the PART clustering analysis
#' @slot sem_sim_org A string indicating the annotation DBI organism to use
#' @slot Gprofiler_results A list of Gprofiler results for the PART clusters
#'
#' @name TimeSeries_Object-class
#' @rdname TimeSeries_Object-class
#' @concept objects
#'
#' @import rstudioapi
#' @importClassesFrom DESeq2 DESeqDataSet
#' @importClassesFrom limma EList
#' @importClassesFrom GOSemSim GOSemSimDATA
#'
#'
#' @exportClass TimeSeries_Object
#'
TimeSeries_Object<-setClass(
  Class='TimeSeries_Object',slots=list(
    count_matrix='list',
    sample_data='data.frame',
    group_names='vector',
    group_colors='vector',
    DE_method='character',
    DE_p_filter='character',
    DE_p_thresh='numeric',
    DE_l2fc_thresh='numeric',
    Gpro_org='character',
    sem_list='GOSemSimDATA',
    DESeq2_obj='DESeqDataSet',
    limma_object='EList',
    DE_results='list',
    PART_l2fc_thresh='numeric',
    PART_results='list',
    sem_sim_org='character',
    Gprofiler_results='list'))


# Class related functions  -------------------
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @title Sample data set-up
#'
#' @description A simple function to read and subset the sample_data created by the user
#'
#'
#' @param path The path that will be read to obtain the sample data
#' @param group_names The group names selected by the user
#'
#' @return the subsetted sample file based on the selected groups
#'
#' @examples
#'  write_example_data_to_dir('PBMC')
#' my_path_data<-'data/PBMC/raw_counts_TS'
#' my_path_sample_dta<-'data/PBMC/sample_file.csv'
#' prep_sample_data(my_path_sample_dta,c('IgM','LPS'))
#'
#' @export
#'
prep_sample_data<-function(path, group_names){
  sample_file<-read.csv(path)
  sample_file<-sample_file[sample_file$group %in% group_names,]

  return(sample_file)

}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @title Prep TS matrix
#'
#' @description Function which prepares and ads the raw count matrix to the timeseries object
#'
#' The function also differs in the preparation method depending onthe type of
#' differential gene expression method selected (DESeq2 or limma) as this affects
#' the input data. Where DESeq2 (RNAseq data) is stored as individual count file while
#' for limma the function expects a csv file with all samples and genes.
#'
#'
#' @param time_object A timeseries object
#' @param path_to_data The path towards the count files.
#' @param limma_from_raw Indicates if raw data is provided for microarrays or if an
#' Elist was given
#' @param limma_id_replace The name of gene id's to use from the 'genes' dataframe
#' in the Elist
#'
#' @return The timeseries object with the raw count matrix added to it
#'
#' @examples
#' write_example_data_to_dir('PBMC')
#' my_path_data<-'data/PBMC/raw_counts_TS'
#' my_path_sample_dta<-'data/PBMC/sample_file.csv'
#'
#' TS_object <- new('TimeSeries_Object',sample_data=prep_sample_data(my_path_sample_dta,c('IgM','LPS')),
#'                  group_names=c('IgM','LPS'),group_colors=c("#e31a1c","#1f78b4"),DE_method='DESeq2',
#'                  DE_p_filter='padj',DE_p_thresh=0.05,DE_l2fc_thresh=1,
#'                  PART_l2fc_thresh=4,sem_sim_org='org.Hs.eg.db',Gpro_org='hsapiens')
#' TS_object <- create_raw_count_matrix(TS_object,my_path_data)
#'
#' @export
create_raw_count_matrix<-function(time_object,path_to_data=NULL,limma_id_replace='GeneName'){

  groups<-slot(time_object,'group_names')
  #Ensures that the order will follow the grouping order
  samp_data<-slot(time_object,'sample_data')
  selected_samples_1<-samp_data$sample[samp_data$group %in% groups[1]]
  selected_samples_2<-samp_data$sample[samp_data$group %in% groups[2]]
  selected_samples<-c(selected_samples_1,selected_samples_2)

  #Prepare the matrix according to the differential expression method (affects input)
  DE_method<-slot(time_object,'DE_method')
  if (DE_method=='DESeq2'){
    final_counts<-prep_RNAseq_matrix(path_to_data,selected_samples)
    slot_storage<-'raw'
  }else if (DE_method=='limma'){
    if(endsWith(path_to_data,'.rds')){
      time_object@limma_object<-readRDS(my_path_data)
    }else{
      time_object@limma_object<-process_microarr_dta_limma(my_path_data)
    }
    Elist<-slot(time_object,'limma_object')
    final_counts<-prep_limma_matrix(Elist_obj=Elist,replace_rows_with = limma_id_replace)
    slot_storage<-'norm' #Stored as norm as the data is expected to be normalized already
  }
  #Re-organize samples to be in order of groups
  final_counts<-final_counts[,selected_samples]


  time_object@count_matrix[[slot_storage]]<-final_counts
  return(time_object)
}

#' @title Read and process microarray data using the limma package
#'
#' @description The function reads in the microarray data, removes the control probes, performs
#' background correction, normalization, and averages over irregular replicate probes (avereps)
#'
#' @param micro_arr_path String indicating the path to a folder containing all and only
#' microrna text
#' @param micro_arr_source The source of the microarray data. Choices can be seen using
#' limma documentation for the function \code{read.maimages}
#' @param green.only Boolean for the green only parameter of \code{read.maimages}
#' @param backg_corr_meth The background correction method with \code{backgroundCorrect}
#' @param back_corr_offset The offset integer to use with \code{backgroundCorrect}
#' @param norm_arrays_meth The normalization method to use with \code{normalizeBetweenArrays}
#' @param ID_used The gene IDs used for the Elist. The name given must match a column in the
#' 'genes' dataframe of the data. If set to NULL, the first column is taken
#'
#' @return The processed and normalized Elist object
#'
#' @importFrom limma read.maimages backgroundCorrect normalizeBetweenArrays avereps
#'
#' @export
process_microarr_dta_limma<-function(micro_arr_path,micro_arr_source='agilent',green.only=TRUE,
                                     backg_corr_meth='normexp',back_corr_offset=16,
                                     norm_arrays_meth='quantile',ID_used='GeneName'){
  #Put file names in vector
  files<-list.files(micro_arr_path)
  raw_files_path <- paste(micro_arr_path, files, sep='/') # add path to filenames

  # load data (see section 4.5 in the limma user guide).
  rawdata <- read.maimages(file=raw_files_path, source=micro_arr_source, green.only=green.only)

  #Rename columns to not include the path to the files
  new_colnames<-gsub(x = colnames(rawdata$E),pattern = paste0(micro_arr_path,'/'),replacement = '')
  colnames(rawdata$E)=new_colnames

  #Remove control
  ctrl_idx <- rawdata@.Data[[4]]$ControlType != 0
  rawdata_noctrl <- rawdata[!ctrl_idx,]

  #Background correction and normalization
  data_bckcor <- backgroundCorrect(rawdata_noctrl, method=backg_corr_meth, offset=back_corr_offset)
  norm_data <- normalizeBetweenArrays(data_bckcor, method=norm_arrays_meth)

  #Create gene ID vector and average the reps
  if(is.null(ID_used)==TRUE){
    ID_vect<-norm_data$genes[,1]
  }else{
    ID_vect<-norm_data$genes[[ID_used]]
  }
  final_dta <- avereps(norm_data, ID=ID_vect)

  return(final_dta)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @title Prepare data matrix from a Elist object
#'
#' @description The function can replace ronames of the matrix by one contained in the 'gene'
#' container of the Elist
#'
#'
#' @param Elist_obj A Elist object
#' @param replace_rows_with Character indicating with column (if any) will be used
#' to replace the rownames of the matrix in the Elist
#'
#' @return count_matrix
#'
#'
#' @export
prep_limma_matrix<-function(Elist_obj,replace_rows_with=NULL){

  count_matrix=Elist_obj$E
  #Replaces rownames with a selection from the genes dataframe of the E list
  if (is.null(replace_rows_with)==FALSE){
    row.names(count_matrix)=Elist_obj$genes[[replace_rows_with]]
  }

  return(count_matrix)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @title prep count matrix form RNAseq count files
#'
#' @description Function requires the path to the individual count files and the samples to
#' be included in the matrix. Files will be read based on the sample they represent
#' the values from the different files are merged into a matrix and returned
#'
#'
#'
#' @param path_to_counts Path to the individual RNAseq count files
#' @param selected_samples The samples to be read
#'
#' @return final_counts The formatted count matrix
#'
#' @examples
#' write_example_data_to_dir('PBMC')
#' my_path_data<-'data/PBMC/raw_counts_TS'
#' my_path_sample_dta<-'data/PBMC/sample_file.csv'
#'
#' TS_object <- new('TimeSeries_Object',sample_data=prep_sample_data(my_path_sample_dta,c('IgM','LPS')),
#'                  group_names=c('IgM','LPS'),group_colors=c("#e31a1c","#1f78b4"),DE_method='DESeq2',
#'                  DE_p_filter='padj',DE_p_thresh=0.05,DE_l2fc_thresh=1,
#'                  PART_l2fc_thresh=4,sem_sim_org='org.Hs.eg.db',Gpro_org='hsapiens')
#' groups<-TS_object@group_names
#' #Ensures that the order will follow the grouping order
#' selected_samples_1<-TS_object@sample_data$sample[TS_object@sample_data$group %in% groups[1]]
#' selected_samples_2<-TS_object@sample_data$sample[TS_object@sample_data$group %in% groups[2]]
#' selected_samples<-c(selected_samples_1,selected_samples_2)
#'
#' #Prepare the matrix according to the differential expression method (affects input)
#' final_counts<-prep_RNAseq_matrix(my_path_data,selected_samples)
#' @export
#'
prep_RNAseq_matrix<-function(path_to_counts,selected_samples){
  final_counts<-data.frame(NULL)
  for(file in list.files(path_to_counts)){
    sample_name<-strsplit(file,'\\.')[[1]][1]
    if (sample_name %in% selected_samples){
      if(endsWith(file,'.csv')==TRUE){
        temp_df<-read.csv(paste0(path_to_counts,'/',file),header=TRUE)
      }else{
        temp_df<-read.delim(paste0(path_to_counts,'/',file),header=FALSE)
      }
      colnames(temp_df)=c('gene_id',sample_name)
      if (nrow(final_counts)==0){
        final_counts<-temp_df
      }else{
        final_counts<-merge(final_counts,temp_df,by='gene_id')
      }
    }
  }
  #Remove any non-genes (starts with underscore)
  final_counts<-final_counts[!startsWith(final_counts$gene_id,'_'),]
  row.names(final_counts)=final_counts$gene_id
  final_counts<-final_counts[,colnames(final_counts)!='gene_id']

  return(final_counts)
}


#' @title Prep example data in directory
#'
#' @description This function writes the example dataset to the 'data' directory
#' This presents the data in a way where users will understand the various files
#' they must provide for the TimeSeriesAnalysis pipeline as well as their format.
#'
#' Having these files written as tab delimited and csv allows users to open them
#' and see the formatting expected by TimeSeriesAnalysis
#'
#' @param example_data Either 'PBMC','MURINE', or 'CELEGANS' to select one of the three example
#' datasets available
#'
#' @return None
#'
#' @examples
#' write_example_data_to_dir('PBMC')
#'
#' @export
write_example_data_to_dir<-function(example_data){
  #If data folder does not exist, create it
  if(!'data' %in% list.files()){
    dir.create('data')
  }
  if(example_data=='PBMC'){
    dir.create('data/PBMC')
    dir.create('data/PBMC/raw_counts_TS')
    for(dta in names(AID_TS_data[['counts']])){
      write.table(AID_TS_data[['counts']][[dta]],paste0('data/PBMC/raw_counts_TS/',dta),quote =FALSE,row.names = FALSE,sep='\t',col.names = FALSE)
    }
    write.csv(AID_TS_data[['sample_dta']],'data/PBMC/sample_file.csv',row.names = FALSE)
  }else if (example_data=='MURINE'){
    dir.create('data/murine')
    dir.create('data/murine/raw_counts_TS')
    for(dta in names(murine_TS_data[['counts']])){
      write.table(murine_TS_data[['counts']][[dta]],paste0('data/murine/raw_counts_TS/',dta),quote =FALSE,row.names = FALSE,sep='\t',col.names = FALSE)
    }
    write.csv(murine_TS_data[['sample_dta']],'data/murine/sample_file.csv',row.names = FALSE)
  }else if(example_data=='CELEGANS'){
    dir.create('data/celegans')
    dir.create('data/celegans/raw_counts_TS')
    for(dta in names(Celegans_TS_data[['counts']])){
      write.table(Celegans_TS_data[['counts']][[dta]],paste0('data/celegans/raw_counts_TS/',dta),quote =FALSE,row.names = FALSE,sep='\t',col.names = FALSE)
    }
    write.csv(Celegans_TS_data[['sample_dta']],'data/celegans/sample_file.csv',row.names = FALSE)
  }else{
    message("Please use 'PBMC','CELEGANS', or 'MURINE' as a parameter to this function.")
  }
}

#' @title Calculate and add semantic similarity to object
#'
#' @description This function calculates the semantic similarity object for an
#' organism. It then stores this within the TimeSeries object and returns it
#'
#' @param object A timeseries object
#' @param ont_sem_sim An ontology for the semantic similarity. EX: 'BP','MF','CC'
#'
#' @return The updated time object
#'
#' @examples
#' TS_object<-create_example_object_for_R()
#' TS_object <- add_semantic_similarity_data(TS_object,ont_sem_sim='BP')
#'
#' @importFrom GOSemSim godata
#'
#' @export
add_semantic_similarity_data<-function(object,ont_sem_sim){
  #Create semantic data
  sem_org<-slot(object,'sem_sim_org')
  object@sem_list <- godata(sem_org, ont=ont_sem_sim, computeIC=TRUE)

  return(object)
}

#' @title Create a subset of count data for R examples (documentation)
#'
#' @description This function creates a subset of the AID count data in order to
#' present quick examples within the documentation
#'
#' @return The smaller count matrix and sample data as a list
#' @examples
#' example_dta<-create_example_data_for_R()
#' @export
create_example_data_for_R<-function(){
  counts_df<-NULL
  for(exp in names(AID_TS_data$counts)){
    if(is.null(counts_df)==TRUE){
      counts_df<-AID_TS_data$counts[[exp]]
    }else{
      counts_df<-cbind(counts_df,AID_TS_data$counts[[exp]][,2])
    }
  }
  row.names(counts_df)=counts_df[,1]
  counts_df<-counts_df[2:ncol(counts_df)]
  colnames(counts_df)=names(AID_TS_data$counts)

  #Subset count dataframe
  counts_df<-counts_df[1:200,]
  return_list<-list(counts=counts_df,sample_data=AID_TS_data$sample_dta)
  return(return_list)
}

#' @title Create a TimeSeries_Object from example data
#'
#' @description This function creates a TimeSeries_Object from example data
#' to use within documentation examples
#'
#' @return The example TimeSeries_Object
#' @examples
#' TS_object<-create_example_object_for_R()
#' @export
create_example_object_for_R<-function(){

  write_example_data_to_dir('PBMC')
  my_path_data<-'data/PBMC/raw_counts_TS'
  my_path_sample_dta<-'data/PBMC/sample_file.csv'

  graph_vect<-c("#e31a1c","#1f78b4")
  names(graph_vect)<-c('IgM','LPS')


  TS_object <- new('TimeSeries_Object',sample_data=prep_sample_data(my_path_sample_dta,c('IgM','LPS')),
                   group_names=c('IgM','LPS'),group_colors=graph_vect,DE_method='DESeq2',
                   DE_p_filter='padj',DE_p_thresh=0.05,DE_l2fc_thresh=1,
                   PART_l2fc_thresh=4,sem_sim_org='org.Hs.eg.db',Gpro_org='hsapiens')
  TS_object <- create_raw_count_matrix(TS_object,my_path_data)
  return(TS_object)
}
