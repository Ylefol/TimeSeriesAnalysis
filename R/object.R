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
#' @slot exp_data A summarized experiment which will contain the sample data and the count matrix
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
#' @slot PART_l2fc_thresh A integer indicating the log(2)foldchange threshold for genes
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
#' @import methods
#' @importClassesFrom DESeq2 DESeqDataSet
#' @importClassesFrom limma EList
#' @importClassesFrom GOSemSim GOSemSimDATA
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#'
#'
#' @exportClass TimeSeries_Object
#'
TimeSeries_Object<-setClass(
  Class='TimeSeries_Object',slots=list(
    exp_data='SummarizedExperiment',
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
#' @param samp_data dataframe containing the sample data
#' @param limma_id_replace The name of gene id's to use from the 'genes' dataframe
#' in the Elist
#'
#' @return The time_object, count matrix, and the slot to save it in
#'
#'
#' @export
create_raw_count_matrix<-function(time_object,samp_data,path_to_data=NULL,limma_id_replace='GeneName'){

  groups<-slot(time_object,'group_names')
  #Ensures that the order will follow the grouping order
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

  return(list(time_object,final_counts,slot_storage))
}

#' @title Add counts and sample data of examples to a TimeSeriesObject
#'
#' @description A function takes an existing TimeSeriesObject and adds the specified
#' example data to the object. Added data is the count matrix and sample data
#'
#'
#' @param time_object A timeseries object
#' @param example_dta Character of either 'PBMC', 'MURINE', or 'CELEGANS'
#' Elist was given
#'
#' @return The timeseries object with the raw count matrix added to it as well as the sample data
#'
#' @examples
#'
#'
#' TS_object <- new('TimeSeries_Object',
#'                  group_names=c('IgM','LPS'),group_colors=c("#e31a1c","#1f78b4"),DE_method='DESeq2',
#'                  DE_p_filter='padj',DE_p_thresh=0.05,DE_l2fc_thresh=1,
#'                  PART_l2fc_thresh=4,sem_sim_org='org.Hs.eg.db',Gpro_org='hsapiens')
#' TS_object <- TS_load_example_data(TS_object,'PBMC')
#'
#' @import SummarizedExperiment
<<<<<<< HEAD
=======
#' @import BiocManager
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import org.Ce.eg.db
#'
>>>>>>> dev
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#'
#' @export
TS_load_example_data<-function(time_object,example_dta){
  if(example_dta=='PBMC'){
    exp_data<-AID_TS_data
  }else if(example_dta=='MURINE'){
    exp_data<-murine_TS_data
  }else if(example_dta=='CELEGANS'){
    exp_data<-Celegans_TS_data
  }else{
    stop("Please use 'PBMC','CELEGANS', or 'MURINE' as a parameter to this function.")
  }
  groups<-slot(time_object,'group_names')
  samp_dta<-data.frame(colData(exp_data))
  count_matrix<-assays(exp_data)$counts
  target_samples<-samp_dta$sample[samp_dta$group %in% groups]

  exp_data<-exp_data[,exp_data$sample %in% target_samples]
  names(assays(exp_data))=c('raw')
  time_object@exp_data<-exp_data

  return(time_object)
}

#' @title Fetches a data matrix from summarized experiments
#'
#' @description Function which extracts a matrix from summarizedexperiments contained
#' within a TimeSeries_object. The matrix is extracted by name.
#'
#'
#' @param time_object A timeseries object
#' @param matrix_name Character for the name of the matrix. 'raw' or 'norm' are the standard
#' names found unless user modifications have been made
#'
#' @return The extracted matrix
#'
#' @examples
#'
#' TS_object <- new('TimeSeries_Object',
#'                  group_names=c('IgM','LPS'),group_colors=c("#e31a1c","#1f78b4"),DE_method='DESeq2',
#'                  DE_p_filter='padj',DE_p_thresh=0.05,DE_l2fc_thresh=1,
#'                  PART_l2fc_thresh=4,sem_sim_org='org.Hs.eg.db',Gpro_org='hsapiens')
#' TS_object <- TS_load_example_data(TS_object,'PBMC')
#' exp_matrix(TS_object,'raw')
#'
#' @importFrom SummarizedExperiment assays
#'
#' @export
exp_matrix<-function(time_object,matrix_name){
  extracted_matrix<-assays(slot(time_object,'exp_data'))[[matrix_name]]
  return(extracted_matrix)
}

#' @title Fetches the sample data
#'
#' @description Function which extracts the sample data from the SummarizedExperiment
#' contained within the TimeSerie_Object. The sample data is extracted and returned in
#' a data frame format
#'
#'
#' @param time_object A timeseries object
#'
#' @return The data.frame format of sample data
#'
#' @examples
#'
#' TS_object <- new('TimeSeries_Object',
#'                  group_names=c('IgM','LPS'),group_colors=c("#e31a1c","#1f78b4"),DE_method='DESeq2',
#'                  DE_p_filter='padj',DE_p_thresh=0.05,DE_l2fc_thresh=1,
#'                  PART_l2fc_thresh=4,sem_sim_org='org.Hs.eg.db',Gpro_org='hsapiens')
#' TS_object <- TS_load_example_data(TS_object,'PBMC')
#' exp_sample_data(TS_object)
#'
#' @importFrom SummarizedExperiment assays
#'
#' @export
<<<<<<< HEAD
exp_sample_data<-function(time_object,matrix_name){
=======
exp_sample_data<-function(time_object){
>>>>>>> dev
  samp_dta<-data.frame(colData(slot(time_object,'exp_data')))
  return(samp_dta)
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


#' @title Adds experiment data in the form of a SummarizedExperiment
#'
#' @description The function creates a SummarizedExperiment from the provided data.
#' It reads the sample data and value files, filters them to contain only the necessary
#' groups and then creates a SummarizedExperiment to then store this into the TimeSeris_Object
#'
#'
#' @param time_object A TimeSeries_Object
#' @param sample_dta_path String which gives the csv path to the sample data
#' @param count_dta_path String which gives the path to the count or microarray data
#' @param limma_ID_replace The name of gene id's to use from the 'genes' dataframe
#' in the Elist
#'
#' @return count_matrix
#'
#' @examples
#' write_example_data_to_dir('PBMC')
#' my_path_data<-'data/PBMC/raw_counts_TS'
#' my_path_sample_dta<-'data/PBMC/sample_file.csv'
<<<<<<< HEAD
=======
#' graph_vect<-c("#e31a1c","#1f78b4")
>>>>>>> dev
#'
#' TS_object <- new('TimeSeries_Object',
#'                  group_names=c('IgM','LPS'),group_colors=graph_vect,DE_method='DESeq2',
#'                  DE_p_filter='padj',DE_p_thresh=0.05,DE_l2fc_thresh=1,
#'                  PART_l2fc_thresh=4,sem_sim_org='org.Hs.eg.db',Gpro_org='hsapiens')
#'
#' TS_object <- add_experiment_data(TS_object,sample_dta_path=my_path_sample_dta,count_dta_path=my_path_data)
#'
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment assays
#'
#' @export
add_experiment_data<-function(time_object,sample_dta_path,count_dta_path,limma_ID_replace='GeneName'){
  groups<-slot(time_object,'group_names')
  sample_data<-prep_sample_data(sample_dta_path,groups)
  sample_data<-sample_data[sample_data$group %in% groups,]
  returned_list<-create_raw_count_matrix(time_object = time_object,path_to_data = count_dta_path,limma_id_replace = limma_ID_replace, samp_data=sample_data)

  time_object<-returned_list[[1]]
  count_matrix<-returned_list[[2]]
  matrix_name<-returned_list[[3]]

  #Add row.names and re-order
  row.names(sample_data)=sample_data$sample
  sample_data<-sample_data[colnames(count_matrix),]

  exp_data<-SummarizedExperiment(assays = count_matrix,colData = sample_data)
  names(assays(exp_data))=matrix_name

  time_object@exp_data<-exp_data

  return(time_object)
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
#' graph_vect<-c("#e31a1c","#1f78b4")
#'
#' TS_object <- new('TimeSeries_Object',
#'                  group_names=c('IgM','LPS'),group_colors=graph_vect,DE_method='DESeq2',
#'                  DE_p_filter='padj',DE_p_thresh=0.05,DE_l2fc_thresh=1,
#'                  PART_l2fc_thresh=4,sem_sim_org='org.Hs.eg.db',Gpro_org='hsapiens')
#'
#' TS_object <- add_experiment_data(TS_object,sample_dta_path=my_path_sample_dta,count_dta_path=my_path_data)
<<<<<<< HEAD
#' groups<-slot(time_object,'group_names')
#' #Ensures that the order will follow the grouping order
#' selected_samples_1<-samp_data$sample[samp_data$group %in% groups[1]]
#' selected_samples_2<-samp_data$sample[samp_data$group %in% groups[2]]
#' selected_samples<-c(selected_samples_1,selected_samples_2)
#'
#' #Prepare the matrix according to the differential expression method (affects input)
#'  final_counts<-prep_RNAseq_matrix(path_to_data,selected_samples)
=======
#' groups<-slot(TS_object,'group_names')
#' #Ensures that the order will follow the grouping order
#' sample_data<-exp_sample_data(TS_object)
#' selected_samples_1<-sample_data$sample[sample_data$group %in% groups[1]]
#' selected_samples_2<-sample_data$sample[sample_data$group %in% groups[2]]
#' selected_samples<-c(selected_samples_1,selected_samples_2)
#'
#' #Prepare the matrix according to the differential expression method (affects input)
#'  final_counts<-prep_RNAseq_matrix(my_path_data,selected_samples)
>>>>>>> dev
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
#'
#' @import SummarizedExperiment
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#'
#' @export
write_example_data_to_dir<-function(example_data){
  #Set object and save names based on requested example data
  if(example_data=='PBMC'){
    full_counts<-assay(AID_TS_data)
    sample_data<-data.frame(AID_TS_data@colData)
    save_name<-'PBMC'
  }else if(example_data=='MURINE'){
    full_counts<-assay(murine_TS_data)
    sample_data<-data.frame(murine_TS_data@colData)
    save_name='murine'
  }else if(example_data=='CELEGANS'){
    full_counts<-assay(Celegans_TS_data)
    sample_data<-data.frame(Celegans_TS_data@colData)
    save_name='celegans'
  }else{
    message("Please use 'PBMC','CELEGANS', or 'MURINE' as a parameter to this function.")
    return(NULL)
  }

  #Create data folder is it does not exist
  if(!'data' %in% list.files()){
    dir.create('data')
  }
  #Create folders to save example data
  dir.create(paste0('data/',save_name))
  dir.create(paste0('data/',save_name,'/raw_counts_TS'))

  #Iterate over data and save it
  for(dta in colnames(full_counts)){
    temp_df<-data.frame(gene_id=row.names(full_counts),samp=full_counts[,dta])
    colnames(temp_df)=c('gene_id',dta)
    write.table(temp_df,paste0('data/',save_name,'/raw_counts_TS/',dta,'.counts'),quote =FALSE,row.names = FALSE,sep='\t',col.names = FALSE)
  }
  write.csv(sample_data,paste0('data/',save_name,'/sample_file.csv'),row.names = FALSE)
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
#'
#' @import SummarizedExperiment
#' @import org.Hs.eg.db
#'
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @export
create_example_data_for_R<-function(){
  counts_df<-NULL
  for(exp in colnames(assays(AID_TS_data)$counts)){
    if(is.null(counts_df)==TRUE){
      counts_df<-assays(AID_TS_data)$counts[,exp]
    }else{
      counts_df<-cbind(counts_df,assays(AID_TS_data)$counts[,exp])
    }
  }
  colnames(counts_df)=colnames(assays(AID_TS_data)$counts)

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
#'
#' @import SummarizedExperiment
#' @import org.Hs.eg.db
#'
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#'
#' @export
create_example_object_for_R<-function(){

  write_example_data_to_dir('PBMC')
  my_path_data<-'data/PBMC/raw_counts_TS'
  my_path_sample_dta<-'data/PBMC/sample_file.csv'

  graph_vect<-c("#e31a1c","#1f78b4")
  names(graph_vect)<-c('IgM','LPS')


  TS_object <- new('TimeSeries_Object',
                   group_names=c('IgM','LPS'),group_colors=graph_vect,DE_method='DESeq2',
                   DE_p_filter='padj',DE_p_thresh=0.05,DE_l2fc_thresh=1,
                   PART_l2fc_thresh=4,sem_sim_org='org.Hs.eg.db',Gpro_org='hsapiens')
  TS_object <- TS_load_example_data(TS_object,'PBMC')
  return(TS_object)
}



