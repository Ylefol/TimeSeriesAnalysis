# library(tibble)
# library(dplyr)

suppressPackageStartupMessages(lapply(packages_for_loading, require, character.only = TRUE))

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
#'
#' @return The timeseries object with the conitional differential expression results
#' added to the DE_results slot of the object.
#'
#' @export
conditional_DE_wrapper<-function(time_object){
  group_names<-time_object@group_names

  if ('conditional' %in% names(time_object@DE_results)){
    message('Conditional differential expression results already exist')
    return(time_object)
  }else{
    time_object@DE_results[['conditional']]<-list()
  }

  all_timepoints<-unique(time_object@sample_data$timepoint)
  for (tp in all_timepoints){
    exp_name<-paste0(group_names[1],'_vs_',group_names[2],'_TP_',tp)
    if (time_object@DE_method=='DESeq2'){
      samps_interest<-time_object@sample_data$sample[time_object@sample_data$timepoint==tp]
      time_object<-DE_using_DESeq2(time_object,group_names,samps_interest,exp_name,main_key='conditional')
    }else if(time_object@DE_method=='limma'){
      time_object<-DE_using_limma(time_object,group_names,exp_name,target_tp=tp,do_temporal=F)
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
#'
#' @return The timeseries object with the temporal differential expression results
#' added to the DE_results slot of the object.
#'
#' @export
#'
temporal_DE_wrapper<-function(time_object){

  if ('temporal' %in% names(time_object@DE_results)){
    message('Temporal differential expression results already exist')
    return(time_object)
  }else{
    time_object@DE_results[['temporal']]<-list()
  }

  all_timepoints<-unique(time_object@sample_data$timepoint)
  for (tp_idx in 1:length(all_timepoints)){
    if ((tp_idx)==length(all_timepoints)){#Stop iteration at second last one
      break
    }
    tps_interest<-c(all_timepoints[tp_idx],all_timepoints[tp_idx+1])
    tp_labels<-paste0('TP_',tps_interest)
    exp_name<-paste0(tp_labels[2],'_vs_',tp_labels[1])

    if (time_object@DE_method=='DESeq2'){
      #Reverse the labels to get the appropriate comparison from the differential gene expression
      rev_tp_labes<-rev(tp_labels)

      samps_interest<-time_object@sample_data$sample[time_object@sample_data$timepoint %in% tps_interest]

      samp_group_1<-time_object@sample_data$sample[time_object@sample_data$timepoint==tps_interest[1]]
      samp_group_2<-time_object@sample_data$sample[time_object@sample_data$timepoint==tps_interest[2]]
      sample_order<-row.names(time_object@DESeq2_obj@colData)
      sample_order<-sample_order[sample_order %in% samps_interest]

      cond_vect<-c()
      for(samp in sample_order){
        if (samp %in% samp_group_1){
          cond_vect<-c(cond_vect,tp_labels[1])
        }else if (samp %in% samp_group_2)
          cond_vect<-c(cond_vect,tp_labels[2])
      }
      cond_vect<-factor(cond_vect,levels = c(tp_labels[2],tp_labels[1]))

      time_object<-DE_using_DESeq2(time_object,rev_tp_labes,sample_order,exp_name,'temporal',condition_factor=cond_vect)
    }else if(time_object@DE_method=='limma'){
      time_object<-DE_using_limma(time_object,group_names=tps_interest,exp_name,do_temporal=T)
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
#' The function then runs through the standard DESeq2 pipeline. It saves the dds
#' object, the raw differential expression results, and the significant differential
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
#' @export
#'
DE_using_DESeq2<-function(time_object,groups,samples_to_use,exp_name,main_key,condition_factor=NULL){

  #Extract significance filtering values
  filter_type<-time_object@DE_p_filter
  pval_thresh<-time_object@DE_p_thresh
  l2fc_thresh<-time_object@DE_l2fc_thresh

  #Subsetting method
  dds<-time_object@DESeq2_obj[,samples_to_use]

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

  res_list<-list(sub_dds=dds,DE_raw_data=DE_results, DE_sig_data=DE_sig)

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
#'
#' @return The timeseries object with the added normalized count matrix
#'
#' @export
#'
normalize_timeSeries_with_deseq2 <- function(time_object){

  groups<-time_object@group_names
  #Assign Condition to all of the columns/samples/patients
  condition <- factor(c(rep(groups[1], nrow(time_object@sample_data[time_object@sample_data$group==groups[1],])),
                        rep(groups[2], nrow(time_object@sample_data[time_object@sample_data$group==groups[2],]))))
  my_matrix<-as.matrix(round(time_object@count_matrix$raw))

  #Create a coldata frame and instantiate the DESeqDataSet
  col_data <- data.frame(row.names=colnames(my_matrix),condition)
  dds <- DESeqDataSetFromMatrix(countData=as.matrix(my_matrix), colData=col_data, design=~condition)
  dds = estimateSizeFactors(dds)

  norm_counts <- as.data.frame(counts(dds,normalized=TRUE))
  norm_counts<- na.omit(norm_counts)

  #Re-orders norm as to have the same order as the raw assay, prevents re-ordering issues
  #Downstream in the pipeline
  # norm_counts<-norm_counts[,colnames(time_object@count_matrix$raw)]
  time_object@count_matrix$norm<-norm_counts
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
#' @export
#'
prep_tp_matrix <- function(time_object,groups_in_ts,target_tp=NULL,do_temporal=FALSE){
  #extract sample data for readability of code
  sample_dta<-time_object@sample_data
  if(do_temporal==FALSE){
    group_ctrl<-sample_dta$sample[sample_dta$group==groups_in_ts[1]]
    group_exp<-sample_dta$sample[sample_dta$group==groups_in_ts[2]]
    tp_samples<-sample_dta$sample[sample_dta$timepoint==target_tp]
    group_1_tp<-group_exp[group_exp %in% tp_samples]
    group_2_tp<-group_ctrl[group_ctrl %in% tp_samples]
  }else{
    group_1_tp<-sample_dta$sample[sample_dta$timepoint==groups_in_ts[2]]
    group_2_tp<-sample_dta$sample[sample_dta$timepoint==groups_in_ts[1]]
  }


  subset_samples<-c(group_1_tp,group_2_tp)

  # tp_limma_subset <- rna_biop_dataa[,colnames(rna_biop_dataa$E) %in% subset_samples]
  microarr_dta<-time_object@limma_object
  g1_vect<-gsub(pattern = F,replacement = 0,x = colnames(microarr_dta$E) %in% group_1_tp)
  g1_vect<-gsub(pattern = T,replacement = 1,x = g1_vect)


  g2_vect<-gsub(pattern = F,replacement = 0,x = colnames(microarr_dta$E) %in% group_2_tp)
  g2_vect<-gsub(pattern = T,replacement = 1,x = g2_vect)
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
#' @export
#'
calculate_EB <- function(micro_arr_dta,matrix_model,comparison){
  if(is.fullrank(matrix_model)==F){
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
DE_using_limma<-function(time_object,group_names,exp_name,target_tp=NULL,do_temporal=F){

  return_list<-prep_tp_matrix(time_object,group_names,target_tp=target_tp,do_temporal=do_temporal)

  my_mm<-return_list[[1]]
  all_samp_IDs<-return_list[[2]]
  my_eb_res<-calculate_EB(time_object@limma_object,my_mm,exp_name)

  if(do_temporal==T){
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
#' @param DE_type Either temporal or conditional for the type of differential gene
#' expression experiment performed
#'
#' @return The updated time object
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

  samples_used_df<-time_object@count_matrix$norm[,samples_used]
  samples_used_df<-samples_used_df[DE_raw$gene_id,]

  DE_raw<-cbind(DE_raw,samples_used_df)

  DE_sig<-DE_raw[abs(DE_raw$log2FoldChange)>time_object@DE_l2fc_thresh,]
  DE_sig<-DE_sig[DE_sig[[time_object@DE_p_filter]]<time_object@DE_p_thresh,]

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

  if(is.null(custom_l2fc_thresh)==T){
    custom_l2fc_thresh<-time_object@PART_l2fc_thresh
  }
  gene_vect<-c()
  for (DE_type in names(time_object@DE_results)){
    DE_res_exps<-time_object@DE_results[[DE_type]]
    for (exp in names(DE_res_exps)){
      genes_above_thresh<-DE_res_exps[[exp]][['DE_sig_data']]
      genes_above_thresh<-genes_above_thresh$gene_id[abs(genes_above_thresh$log2FoldChange)>=custom_l2fc_thresh]
      gene_vect<-c(gene_vect,genes_above_thresh)
    }
  }
  gene_vect<-unique(gene_vect)
  return(gene_vect)
}
