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
#' @export
#'
find_most_variable_cluster<-function(time_object,mean_ts_data){

  #Retrieve clusters which have gprofiler results
  clusts_interest<-names(time_object@Gprofiler_results)
  mean_ts_data<-mean_ts_data[mean_ts_data$cluster %in% clusts_interest,]

  dif_clust<-c()
  for(clust in unique(mean_ts_data$cluster)){
    dif_vect<-c()
    for(tp in unique(mean_ts_data$timepoint)){
      val_exp<-mean_ts_data$trans_mean[mean_ts_data$timepoint==tp & mean_ts_data$cluster==clust & mean_ts_data$group==TS_object@group_names[1]]
      val_control<-mean_ts_data$trans_mean[mean_ts_data$timepoint==tp & mean_ts_data$cluster==clust & mean_ts_data$group==TS_object@group_names[2]]

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
#' @export
#'
create_DEG_df<-function(time_object){
  #Get DEGs
  DEG_amount<-c()
  names_exp<-c()
  gene_vect<-c()
  for(DE_type in names(time_object@DE_results)){
    for(exp in names(time_object@DE_results[[DE_type]])){
      gene_vect<-c(gene_vect,time_object@DE_results[[DE_type]][[exp]][['DE_sig_data']]$gene_id)
      DEG_amount<-c(DEG_amount,nrow(time_object@DE_results[[DE_type]][[exp]][['DE_sig_data']]))
      exp<-paste0(exp,' (',DE_type,')')
      names_exp<-c(names_exp,exp)
    }
  }
  DEG_amount<-c(DEG_amount,paste0('**',length(unique(gene_vect)),'**'))
  names_exp<-c(names_exp,'**total unique genes**')
  DEG_summary<-data.frame(`experiment name`=names_exp,`number of DEGs`=DEG_amount)

  return(DEG_summary)
}
