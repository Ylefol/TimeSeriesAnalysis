# library(tibble)
# library(dplyr)

# suppressPackageStartupMessages(lapply(packages_for_loading, require, character.only = TRUE))

# wrapper functions  -------------------

#' @title Create conditional DE results
#'
#' @description A wrapper function which performs conditional differential gene expression
#' analysis for every timepoint in the timeseries object
#'
#' Conditional differential gene expression experiments involves the comparison
#' of the experiment vs the control. This is done at every timepoint.
#'
#' The function calls the necessary sub functions for the differential expression
#' analysis based on the type needed (DESeq2 for RNAseq or limma for microarray)
#'
#'
#' @param time_object A timeseries object containing a DESeq2_obj
#' @param vignette_run Boolean indicating if a run is for vignettes or not.
#'
#' @return The timeseries object with the conitional differential expression results
#' added to the DE_results slot of the object.
#'
#' @examples
#' TS_object<-create_example_object_for_R()
#' TS_object<-normalize_timeSeries_with_deseq2(time_object=TS_object)
#' TS_object<-conditional_DE_wrapper(TS_object,vignette_run=TRUE)
#'
#' @export
conditional_DE_wrapper<-function(time_object,vignette_run=FALSE){

  if(vignette_run==TRUE){
    time_object@DE_results$conditional<-PBMC_pre_loaded$DE_results$conditional
    return(time_object)
  }

  group_names<-slot(time_object,'group_names')
  DE_res<-slot(time_object,'DE_results')
  DE_meth<-slot(time_object,'DE_method')
  samp_dta_full<-exp_sample_data(time_object)
  if ('conditional' %in% names(DE_res)){
    message('Conditional differential expression results already exist')
    return(time_object)
  }else{
    DE_res[['conditional']]<-list()
  }

  all_timepoints<-unique(samp_dta_full$timepoint)
  for (tp in all_timepoints){
    exp_name<-paste0(group_names[1],'_vs_',group_names[2],'_TP_',tp)
    if (DE_meth=='DESeq2'){
      samps_interest<-samp_dta_full$sample[samp_dta_full$timepoint==tp]
      time_object<-DE_using_DESeq2(time_object,group_names,samps_interest,exp_name,main_key='conditional')
    }else if(DE_meth=='limma'){
      time_object<-DE_using_limma(time_object,group_names,exp_name,target_tp=tp,do_temporal=FALSE)
    }


  }
  return(time_object)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @title Create temporal DE results
#'
#' @description Wrapper function which performs the temporal differential gene expression analyses
#'
#' The temporal differential analyses is the comparison of each timepoint irrelevant
#' of their conition (experiment or control). Timepoints are compared with the next
#' immediate timepoint. For example, if there are three timepoints, two temporal
#' differential gene expressions will be performed: TP2 vs TP1 and TP3 vs TP2.
#'
#' The function subsets the DESeq2_obj based on the necessary samples and adjusts
#' the condition of the object to reflect the comparison being done
#'
#'
#' @param time_object A timeseries object containing a DESeq2_obj
#' @param do_all_combinations Allows for all temporal combinations to be done instead
#' of just sequential comparison. ex: do TP2vsTP1, TP3vsTP2, AND TP3vsTP1. In a normal instance
#' only the first two comparison of the example would be run.
#' @param vignette_run Boolean indicating if a run is for vignettes or not.
#'
#' @return The timeseries object with the temporal differential expression results
#' added to the DE_results slot of the object.
#'
#' @examples
#' TS_object<-create_example_object_for_R()
#' TS_object<-normalize_timeSeries_with_deseq2(time_object=TS_object)
#' TS_object<-temporal_DE_wrapper(TS_object,do_all_combinations=TRUE,vignette_run=TRUE)
#'
#' @export
#'
temporal_DE_wrapper<-function(time_object,do_all_combinations=FALSE,vignette_run=FALSE){

  if(vignette_run==TRUE){
    time_object@DE_results$temporal<-PBMC_pre_loaded$DE_results$temporal
    return(time_object)
  }

  DE_res<-slot(time_object,'DE_results')
  samp_dta_full<-exp_sample_data(time_object)

  if ('temporal' %in% names(DE_res)){
    message('Temporal differential expression results already exist')
    return(time_object)
  }else{
    DE_res[['temporal']]<-list()
  }
  #Sort the timepoints to ensure that they are in sequential order
  all_timepoints<-sort(unique(samp_dta_full$timepoint))
  if(do_all_combinations==TRUE){
    possible_TPs<-expand.grid(all_timepoints,all_timepoints)
    for (row in 1:nrow(possible_TPs)){
      row<-as.character(row)
      if(possible_TPs[row,'Var1'] <= possible_TPs[row,'Var2']){
        possible_TPs<-possible_TPs[!row.names(possible_TPs) %in% c(row),]
      }
    }
  }else{ #Do all combinations is False
    possible_TPs<-data.frame(NULL)
    for (tp_idx in 1:length(all_timepoints)){
      if ((tp_idx)==length(all_timepoints)){#Stop iteration at second last one
        break
      }
      found_tps<-rev(c(all_timepoints[tp_idx],all_timepoints[tp_idx+1]))
      if(nrow(possible_TPs)==0){
        possible_TPs<-data.frame('Var1'=found_tps[1],'Var2'=found_tps[2])
      }else{
        possible_TPs<-rbind(possible_TPs,found_tps)
      }
    }
  }
  DE_meth<-slot(time_object,'DE_method')
  #Iterate over possible comparison and run them
  for(exp in 1:nrow(possible_TPs)){
    tps_interest<-as.numeric(possible_TPs[exp,])
    tp_labels<-paste0('TP_',tps_interest)
    exp_name<-paste0(tp_labels[1],'_vs_',tp_labels[2])
    if (DE_meth=='DESeq2'){

      samps_interest<-samp_dta_full$sample[samp_dta_full$timepoint %in% tps_interest]
      samp_group_1<-samp_dta_full$sample[samp_dta_full$timepoint==tps_interest[2]]
      samp_group_2<-samp_dta_full$sample[samp_dta_full$timepoint==tps_interest[1]]
      DESeq2_obj<-slot(time_object,'DESeq2_obj')
      DESeq2_colData<-slot(DESeq2_obj,'colData')
      sample_order<-row.names(DESeq2_colData)
      sample_order<-sample_order[sample_order %in% samps_interest]

      #The labels show the highest tp first then the lowest. To build the comparison
      #Vectors we need to put the lowest as sample 1 and highest as sample 2
      cond_vect<-c()
      for(samp in sample_order){
        if (samp %in% samp_group_1){
          cond_vect<-c(cond_vect,tp_labels[2])
        }else if (samp %in% samp_group_2)
          cond_vect<-c(cond_vect,tp_labels[1])
      }
      cond_vect<-factor(cond_vect,levels = c(tp_labels[1],tp_labels[2]))
      time_object<-DE_using_DESeq2(time_object,tp_labels,sample_order,exp_name,'temporal',condition_factor=cond_vect)

    }else if(DE_meth=='limma'){
      time_object<-DE_using_limma(time_object,group_names=tps_interest,exp_name,do_temporal=TRUE)
    }

  }

  return(time_object)
}


# DESeq2 functions  -------------------

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @title DE with DESeq2
#'
#' @description Function which performs the differential gene expression experiment using
#' DESeq2
#'
#' The function subsets the DESeq2_obj if necessary and adjust the condition in the
#' event of a temporal differential gene expression analysis.
#'
#' The function then runs through the standard DESeq2 pipeline. It saves the raw
#' differential expression results, and the significant differential
#' expression results.
#'
#'
#' @param time_object A timeseries object
#' @param groups The groups being used (varies from conditional to temporal)
#' @param samples_to_use Samples to use in the comparisons
#' @param exp_name The name of the experiment, it will be use to store the results
#' @param main_key Either conditional or temporal to indicate which type of
#' differential gene expression experiment is being performed
#' @param condition_factor A factor containing the new condition used for the
#' differeital gene expression analysis
#'
#' @return The timeseries object updated with the results from the experiment
#'
#' @examples
#' TS_object<-create_example_object_for_R()
#' TS_object <- normalize_timeSeries_with_deseq2(time_object=TS_object)
#'
#' #DE for a single timepoint
#' group_names<-slot(TS_object,'group_names')
#' tp<-'1'
#' exp_name<-paste0(group_names[1],'_vs_',group_names[2],'_TP_',tp)
#' sample_data<-exp_sample_data(TS_object)
#' samps_interest<-sample_data$sample[sample_data$timepoint==tp]
#' TS_object<-DE_using_DESeq2(TS_object,group_names,samps_interest,exp_name,main_key='conditional')
#'
#' @importFrom DESeq2 DESeq results
#' @import BiocGenerics
#'
#' @export
#'
DE_using_DESeq2<-function(time_object,groups,samples_to_use,exp_name,main_key,condition_factor=NULL){

  #Extract significance filtering values
  filter_type<-slot(time_object,'DE_p_filter')
  pval_thresh<-slot(time_object,'DE_p_thresh')
  l2fc_thresh<-slot(time_object,'DE_l2fc_thresh')

  #Subsetting method
  dds<-slot(time_object,'DESeq2_obj')[,samples_to_use]

  if(main_key=='temporal'){
    dds@colData@listData[['condition']]<-condition_factor
  }

  dds <- DESeq(dds)

  # plot_deseq2_PCAs(dds,folder_n)
  #Change the format into something that can be used by R (dds --> res)
  DE_results <- results(dds, contrast = c("condition", groups[1], groups[2]))

  #Orders the results via adjusted p-value
  DE_results <- DE_results[order(DE_results['padj']), ]

  #Merges the DE data with the patient data, will be saved as DE_raw_data
  DE_results <- merge(as.data.frame(DE_results), as.data.frame(counts(dds,normalized=TRUE)), by= "row.names", sort=FALSE)

  names(DE_results)[1] <- "gene_id"
  DE_results<- na.omit(DE_results)

  DE_sig<-DE_results[DE_results[filter_type]<pval_thresh & abs(DE_results$log2FoldChange)>l2fc_thresh,]

  # res_list<-list(sub_dds=dds,DE_raw_data=DE_results, DE_sig_data=DE_sig)
  res_list<-list(DE_raw_data=DE_results, DE_sig_data=DE_sig)

  time_object@DE_results[[main_key]][[exp_name]]<-res_list


  return(time_object)

}


#' @title Normalize with DESeq2
#'
#' @description Function which normalizes the raw count matrix from a timeseries object
#' using the DESeq2 normalization.
#'
#' @param time_object The time object with the normalized count matrix included,
#' as well as the DESeq2 object
#' @param controlGenes The controlGenes parameter to be passed to DESeq's
#' estimateSizeFactors function. It allows to normalize on specific identifiers.
#'
#' @return The timeseries object with the added normalized count matrix
#'
#' @examples
#' TS_object<-create_example_object_for_R()
#' TS_object <- normalize_timeSeries_with_deseq2(time_object=TS_object)
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors counts
#' @importFrom SummarizedExperiment assays
#'
#' @export
#'
normalize_timeSeries_with_deseq2 <- function(time_object,controlGenes=NULL){

  groups<-slot(time_object,'group_names')
  samp_dta_full<-exp_sample_data(time_object)
  count_matrix<-exp_matrix(time_object,'raw')

  #Set levels to the inverse of the submitted group names, therefore control is the reference
  condition<-factor(samp_dta_full$group,levels=rev(time_object@group_names))

  my_matrix<-as.matrix(round(count_matrix))

  #Create a coldata frame and instantiate the DESeqDataSet
  col_data <- data.frame(row.names=colnames(my_matrix),condition)
  dds <- DESeqDataSetFromMatrix(countData=as.matrix(my_matrix), colData=col_data, design=~condition)
  if(is.null(controlGenes)==TRUE){
    dds = estimateSizeFactors(object=dds)
  }else{
    numerical_control_index=which(rownames(as.matrix(my_matrix)) %in% controlGenes)
    dds = estimateSizeFactors(object=dds,controlGenes=numerical_control_index)#Ability to input custom sequences
  }



  norm_counts <- as.data.frame(counts(dds,normalized=TRUE))
  norm_counts<- na.omit(norm_counts)

  assays(time_object@exp_data)$norm<-norm_counts
  time_object@DESeq2_obj<-dds
  return(time_object)
}


# limma functions  -------------------

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


#' @title limma timepoint matrix
#'
#' @description Function which prepares a timepoint model matrix for limma datasets
#'
#' The function creates either a conditional matrix (do_temporal==FALSE) where
#' the model matrix will compare both conditions ex: WT vs Knockout
#' Or it will create a model for a temporal comparison, where it will compare the
#' two timpoints contained in 'groups_in_ts'
#'
#' The function is meant to be used within a wrapper function
#'
#' @param time_object A time series object
#' @param groups_in_ts Vector of the two groups/conditionals
#' @param target_tp Vector containing the two timepoints to compare
#' @param do_temporal boolean indication if the model matrix should be temporal or conditional
#'
#' @return a list contianing the model matrix and the subsets used
#'
#' @importFrom stats model.matrix
#'
#' @export
#'
prep_tp_matrix <- function(time_object,groups_in_ts,target_tp=NULL,do_temporal=FALSE){
  #extract sample data for readability of code
  samp_dta_full<-exp_sample_data(time_object)
  if(do_temporal==FALSE){
    group_ctrl<-samp_dta_full$sample[samp_dta_full$group==groups_in_ts[1]]
    group_exp<-samp_dta_full$sample[samp_dta_full$group==groups_in_ts[2]]
    tp_samples<-samp_dta_full$sample[samp_dta_full$timepoint==target_tp]
    group_1_tp<-group_exp[group_exp %in% tp_samples]
    group_2_tp<-group_ctrl[group_ctrl %in% tp_samples]
  }else{
    group_1_tp<-samp_dta_full$sample[samp_dta_full$timepoint==groups_in_ts[2]]
    group_2_tp<-samp_dta_full$sample[samp_dta_full$timepoint==groups_in_ts[1]]
  }


  subset_samples<-c(group_1_tp,group_2_tp)

  # tp_limma_subset <- rna_biop_dataa[,colnames(rna_biop_dataa$E) %in% subset_samples]
  microarr_dta<-slot(time_object,'limma_object')
  g1_vect<-gsub(pattern = FALSE,replacement = 0,x = colnames(microarr_dta$E) %in% group_1_tp)
  g1_vect<-gsub(pattern = TRUE,replacement = 1,x = g1_vect)


  g2_vect<-gsub(pattern = FALSE,replacement = 0,x = colnames(microarr_dta$E) %in% group_2_tp)
  g2_vect<-gsub(pattern = TRUE,replacement = 1,x = g2_vect)
  matrix_df<-data.frame(Group1=g1_vect, Group2=g2_vect)
  mm <- model.matrix(~Group1+Group2,matrix_df)
  colnames(mm)=c('X.Intercept','G1','G2')

  return(list(mm,subset_samples))
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @title Empirical bayesan analysis
#'
#' @description Function which uses a model matrix to calculate the empirical Bayesian results
#' from a provided microarray normalized and corrected Elist.
#'
#' @param micro_arr_dta An Elist object containing the normalized and corrected microarray data
#' @param matrix_model the matrix model as created by the \code{prep_matrix()}
#' @param comparison The name of the comparison being run
#'
#' @return The empirical Bayesian results
#'
#' @importFrom limma is.fullrank lmFit makeContrasts contrasts.fit eBayes
#'
#' @export
#'
calculate_EB <- function(micro_arr_dta,matrix_model,comparison){
  if(is.fullrank(matrix_model)==FALSE){
    matrix_model<-matrix_model[,1:2]
    #Perform linear regression and set-up contrast
    limma_res <- lmFit(micro_arr_dta, matrix_model)
    ctrst <- makeContrasts(G2, levels = colnames(coef(limma_res)))
  }else{
    #Perform linear regression and set-up contrast
    limma_res <- lmFit(micro_arr_dta, matrix_model)
    ctrst <- makeContrasts(G2 - G1, levels = colnames(coef(limma_res)))
  }

  colnames(ctrst) <- c(comparison)
  #Compute Contrasts From Linear Model Fit
  ctrst_res <- contrasts.fit(limma_res, ctrst)
  #Empirical Bayes Statistics For Differential Expression
  eb_res <- eBayes(ctrst_res)
  return (eb_res)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @title DE using limma
#'
#' @description Function which performs differential gene expression for microarray data
#' using limma
#'
#' This function serves as a wrapper function which calls the functions to create
#' the model matrix, the empirical bayesisan calculation (EB), and the function to convert
#' the EB results to dataframe format for the raw and significant differential expression
#' data
#'
#' @param time_object A time series object
#' @param group_names The group names to be used (either conditional or temporal names)
#' @param exp_name The name of the expriment that will be used to store the results
#' @param target_tp The timepoint being targeted in the case of conditional analysis
#' @param do_temporal If the analysis is temporal or conditional
#'
#' @return The updated timeseries object
#'
#' @export
#'
DE_using_limma<-function(time_object,group_names,exp_name,target_tp=NULL,do_temporal=FALSE){

  return_list<-prep_tp_matrix(time_object,group_names,target_tp=target_tp,do_temporal=do_temporal)

  my_mm<-return_list[[1]]
  all_samp_IDs<-return_list[[2]]
  my_eb_res<-calculate_EB(slot(time_object,'limma_object'),my_mm,exp_name)

  if(do_temporal==TRUE){
    DE_type<-'temporal'
  }else{
    DE_type<-'conditional'
  }
  time_object<-convert_eb_res_to_DE_results(time_object,my_eb_res,all_samp_IDs,exp_name,DE_type)

  return(time_object)
}



### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @title Format limma DE results
#'
#' @description Function which converts the empirical bayesian results into data frame format
#' for storage in a DE_results object.
#'
#' The function also formats the data in regards to the column names in order
#' to ensure downstream compatibility.
#'
#' @param time_object A timeseries object
#' @param eb_res The empirical bayesian results
#' @param samples_used The samples used to obtain the eb_res object
#' @param exp_name The name which will be given to the results obtained by the function.
#' @param DE_type Either temporal or conditional for the type of differential gene
#' expression experiment performed
#'
#' @return The updated time object
#'
#' @importFrom limma topTable
#' @importFrom dplyr rename "%>%"
#'
#' @export
#'
convert_eb_res_to_DE_results<-function(time_object,eb_res,samples_used,exp_name,DE_type){
  eb_res %>%
    topTable(coef=1, number = Inf)  -> DE_raw
  DE_raw %>%
    dplyr::rename(
      gene_id = GeneName,
    ) -> DE_raw

  #Filter for columns that we want
  DE_raw <- DE_raw[,c('gene_id','logFC','AveExpr','P.Value','adj.P.Val')]
  #Change column names to match DE results that have been created before to allow
  #these results to be sent to the rest of the pipeline
  colnames(DE_raw)=c('gene_id','log2FoldChange','baseMean','pvalue','padj')

  DE_raw =  DE_raw[order(DE_raw$gene_id,DE_raw$padj),]
  DE_raw = DE_raw[ !duplicated(DE_raw$gene_id), ]

  samples_used_df<-exp_matrix(time_object,'norm')[,samples_used]
  samples_used_df<-samples_used_df[DE_raw$gene_id,]

  DE_raw<-cbind(DE_raw,samples_used_df)

  DE_sig<-DE_raw[abs(DE_raw$log2FoldChange)>slot(time_object,'DE_l2fc_thresh'),]
  DE_sig<-DE_sig[DE_sig[[slot(time_object,'DE_p_filter')]]<slot(time_object,'DE_p_thresh'),]

  res_list<-list(sub_eb=eb_res,DE_raw_data=DE_raw, DE_sig_data=DE_sig)

  if (DE_type=='conditional'){
    time_object@DE_results$conditional[[exp_name]]<-res_list
  }else if (DE_type=='temporal'){
    time_object@DE_results$temporal[[exp_name]]<-res_list
  }

  return(time_object)
}



# General functions  -------------------
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


#' @title Select genes based on L2FC/FC
#'
#' @description Function which goes through all differential gene expression experiments, both
#' conditional and temporal, and extracts all significant genes that have a
#' absolute log2FoldChange greater than the custom_l2fc_thresh
#'
#'
#' @param time_object A timeseries object
#' @param custom_l2fc_thresh A value indicating the log2FoldChange threshold, can be NULL
#'
#' @return A vector of significant genes whose absolute log2FoldChange is greater
#' than the inputed custom threhsold.
#'
#' @export
#'
select_genes_with_l2fc<-function(time_object,custom_l2fc_thresh=NULL){
  DE_res<-slot(time_object,'DE_results')
  if(is.null(custom_l2fc_thresh)==TRUE){
    custom_l2fc_thresh<-slot(time_object,'PART_l2fc_thresh')
  }
  gene_vect<-c()
  for (DE_type in names(DE_res)){
    DE_res_exps<-DE_res[[DE_type]]
    for (exp in names(DE_res_exps)){
      genes_above_thresh<-DE_res_exps[[exp]][['DE_sig_data']]
      genes_above_thresh<-genes_above_thresh$gene_id[abs(genes_above_thresh$log2FoldChange)>=custom_l2fc_thresh]
      gene_vect<-c(gene_vect,genes_above_thresh)
    }
  }
  gene_vect<-unique(gene_vect)
  return(gene_vect)
}
