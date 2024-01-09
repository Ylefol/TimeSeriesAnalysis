#' @title Identify cluster with the most variability
#'
#' @description Identifies the most variable cluster based on the mean trajectory
#' of each cluster. This cluster is identified for the purpose of illustration
#' in an rmarkdown file.
#'
#' @param time_object A timeseries object
#' @param mean_ts_data A dataframe containing the mean trajectory for each cluster
#'
#'
#'
#' @return The name of the most variable cluster
#'
#' @examples
#' TS_object<-create_example_object_for_R()
#' TS_object <- normalize_timeSeries_with_deseq2(time_object=TS_object)
#' #Perform conditional differential gene expression analysis
#' TS_object<-conditional_DE_wrapper(TS_object,vignette_run=TRUE)
#' TS_object<-temporal_DE_wrapper(TS_object,do_all_combinations=TRUE,vignette_run=TRUE)
#' #Extract genes for PART clustering based on defined log(2)foldChange threshold
#' signi_genes<-select_genes_with_l2fc(TS_object)
#'
#' #Use all samples, but implement a custom order. In this case it is reversed
#' sample_data<-exp_sample_data(TS_object)
#' TS_groups<-slot(TS_object,'group_names')
#' samps_2<-sample_data$sample[sample_data$group==TS_groups[2]]
#' samps_1<-sample_data$sample[sample_data$group==TS_groups[1]]
#'
#' #Create the matrix that will be used for PART clustering
#' TS_object<-prep_counts_for_PART(object=TS_object,target_genes=signi_genes,scale=TRUE,target_samples=c(samps_2,samps_1))
#' TS_object<-compute_PART(TS_object,part_recursion=10,part_min_clust=10,dist_param="euclidean", hclust_param="average",vignette_run=TRUE)
#' TS_object<-run_gprofiler_PART_clusters(TS_object,vignette_run=TRUE) #Run the gprofiler analysis
#' ts_data<-calculate_cluster_traj_data(TS_object,scale_feat=TRUE) #Calculate scaled gene values for genes of clusters
#' mean_ts_data<-calculate_mean_cluster_traj(ts_data) #Calculate the mean scaled values for each cluster
#' target_clust<-find_most_variable_cluster(TS_object,mean_ts_data)
#'
#' @export
#'
find_most_variable_cluster<-function(time_object,mean_ts_data){

  #Retrieve clusters which have gprofiler results
  clusts_interest<-names(slot(time_object,'Gprofiler_results'))
  if(is.null(clusts_interest)==TRUE){
    target_clust<-unique(time_object@PART_results$cluster_map$cluster)[1]
    return(target_clust)
  }
  mean_ts_data<-mean_ts_data[mean_ts_data$cluster %in% clusts_interest,]

  dif_clust<-c()
  for(clust in unique(mean_ts_data$cluster)){
    dif_vect<-c()
    for(tp in unique(mean_ts_data$timepoint)){
      val_exp<-mean_ts_data$trans_mean[mean_ts_data$timepoint==tp & mean_ts_data$cluster==clust & mean_ts_data$group==time_object@group_names[1]]
      val_control<-mean_ts_data$trans_mean[mean_ts_data$timepoint==tp & mean_ts_data$cluster==clust & mean_ts_data$group==time_object@group_names[2]]

      dif_val<-abs(val_exp-val_control)
      dif_vect<-c(dif_vect,dif_val)
    }
    dif_clust_val<-sum(dif_vect)
    dif_clust<-c(dif_clust,dif_clust_val)
  }
  names(dif_clust)=unique(mean_ts_data$cluster)
  target_clust<-names(dif_clust)[unname(dif_clust)==max(dif_clust)]
  target_clust<-target_clust[1]
  return(target_clust)

}


#' @title Create summary table for DEGs
#'
#' @description Creates a table which summariezes the conditional and temporal
#' differetial gene expression results. The table shows the number of significant
#' DEGs found in each differential gene experiment performed. The table is created
#' for a Rmarkdown format report
#'
#' @param  time_object A timeseries object
#'
#' @return The DEG table summarizing the results
#'
#' @examples
#' TS_object<-create_example_object_for_R()
#' TS_object <- normalize_timeSeries_with_deseq2(time_object=TS_object)
#' #Perform conditional differential gene expression analysis
#' TS_object<-conditional_DE_wrapper(TS_object,vignette_run=TRUE)
#' TS_object<-temporal_DE_wrapper(TS_object,do_all_combinations=FALSE,vignette_run=TRUE)
#' DEG_df<-create_DEG_df(TS_object)
#'
#' @export
#'
create_DEG_df<-function(time_object){
  DE_res<-slot(time_object,'DE_results')

  DEG_amount<-c()
  names_exp<-c()
  gene_vect<-c()
  for(DE_type in names(DE_res)){
    for(exp in names(DE_res[[DE_type]])){
      gene_vect<-c(gene_vect,DE_res[[DE_type]][[exp]][['DE_sig_data']]$gene_id)
      DEG_amount<-c(DEG_amount,nrow(DE_res[[DE_type]][[exp]][['DE_sig_data']]))
      exp<-paste0(exp,' (',DE_type,')')
      names_exp<-c(names_exp,exp)
    }
  }
  DEG_amount<-c(DEG_amount,paste0('**',length(unique(gene_vect)),'**'))
  names_exp<-c(names_exp,'**total unique genes**')
  DEG_summary<-data.frame(`experiment name`=names_exp,`number of DEGs`=DEG_amount)

  return(DEG_summary)
}
