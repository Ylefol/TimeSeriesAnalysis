

#' @title Read and subset sample data
#'
#' @description Function which reads in the sample data and subsets it for the
#' group names provided
#'
#' @param samp_dta_path The path to the csv of the sample data
#' @param group_names A vector of group_names that will be used to subset the sample
#' data file.
#'
#' @returns The sample_data
#'
#' @export
set_up_sample_dta<-function(samp_dta_path,group_names){

  sample_data<-read.csv(samp_dta_path)
  #subset sample_data for groups of interest
  sample_data<-sample_data[sample_data$group %in% group_names,]
  #sort for the subsequent normalization
  new_sample_data<-data.frame()
  for(group in group_names){
    if(nrow(new_sample_data)==0){
      new_sample_data<-sample_data[sample_data$group==group,]
    }else{
      new_sample_data<-rbind(new_sample_data,sample_data[sample_data$group==group,])
    }
  }
  sample_data<-new_sample_data

  return(sample_data)
}




#' @title Read and process count data
#'
#' @description This function will read in the necessary counts files based on the
#' sample_data file provided. It will look for (in the path to counts) for files
#' which are named after the sample names in the sample data. The count matrix will
#' then be used to create a deseq2 object which will then be normalized.
#'
#' @param path_to_counts The path to the csv of the counts
#' @param sample_data The sample data
#'
#' @import DESeq2
#'
#' @returns A list containing the deseq2 object and the normalized counts matrix
#'
#' @export
create_matrix_and_dds<-function(path_to_counts,sample_data){
  my_counts<-prep_RNAseq_matrix(path_to_counts,sample_data$sample)
  my_counts<-my_counts[,sample_data$sample]
  condition<-factor(sample_data$group,levels=group_names)
  col_data <- data.frame(row.names=colnames(my_counts),condition)
  dds <- DESeqDataSetFromMatrix(countData=as.matrix(my_counts), colData=col_data, design=~condition)
  dds = estimateSizeFactors(object=dds)
  norm_counts <- as.data.frame(counts(dds,normalized=TRUE))
  return(list(dds,norm_counts))
}


#' @title Extract significant genes from several TimeSeriesObjects
#'
#' @description Function which loads given TimeSeriesObjects and extracts the
#' significant genes which have an absolute log 2 fold change greater than the
#' provided threshold (l2fc_thresh). All significant genes are appended to a single
#' vector, duplicates are removed from the vector.
#'
#' @param TS_object_paths A list of path to TimeSeriesObjects of interest
#' @param l2fc_thresh A numerical values for the log2foldchange threshold
#'
#'
#' @returns A vector of the significant genes found
#'
#' @export
get_significant_genes<-function(TS_object_paths,l2fc_thresh){
  target_genes<-c()
  for(name in names(TS_object_paths)){
    load(TS_object_paths[[name]])
    if(length(target_genes)==0){
      target_genes<-select_genes_with_l2fc(TS_object,l2fc_thresh)
    }else{
      target_genes<-c(target_genes,select_genes_with_l2fc(TS_object,l2fc_thresh))
    }
    target_genes<-unique(target_genes)
  }
  return(target_genes)
}


#' @title Creates trajectory plots for clusters of interest
#'
#' @description Function which creates trajectory plots of declared clusters of interest.
#' The function takes in the object, calculates the trajectory data and mean trajectory data
#' for all clusters. It then subsets these values to only contain the clusters of interest.
#' It then submits this information to a plotting function which will create
#' a plot with the desired clusters as well as respecting the other parameters given.
#'
#' @param object A TimeSeriesObject
#' @param target_clusters A vector of clusters of interest
#' @param num_clusters_per_row An integer stating how many clusters will be per row.
#' Note that the number of panels for 'a cluster' is equal to the number of groups in
#' the object.
#' @param log_TP A boolean which states if the timepoints should be log-transformed
#' @param scale If the data should be scaled for the trajectories. By default TRUE.
#'
#' @returns The ggplot2 object
create_cluster_traj_plot<-function(object,target_clusters,num_clusters_per_row,log_TP,scale=TRUE){

  ts_data<-calculate_cluster_traj_data(object,scale_feat=scale) #Calculate scaled gene values for genes of clusters
  if(scale==TRUE){
    yaxis_name='scaled mean counts'
  }else{
    yaxis_name='mean counts'
  }
  mean_ts_data<-calculate_mean_cluster_traj(ts_data) #Calculate the mean scaled values for each cluster

  sub_ts_data<-ts_data[ts_data$cluster%in%target_clusters,]
  sub_ts_data<-sub_ts_data[order(match(sub_ts_data$cluster,target_clusters)),]

  sub_ts_means<-mean_ts_data[mean_ts_data$cluster%in%target_clusters,]
  sub_ts_means<-sub_ts_means[order(match(sub_ts_means$cluster,target_clusters)),]

  sub_ts_data$labels<-factor(sub_ts_data$labels,levels = unique(sub_ts_data$labels))
  sub_ts_means$labels<-factor(sub_ts_means$labels,levels = unique(sub_ts_means$labels))

  cluster_num<-length(target_clusters)
  number_rows<-ceiling(cluster_num/num_clusters_per_row)

  if(log_TP==TRUE){
    sub_ts_data$log10_timepoint<-log10(sub_ts_data$timepoint)
    sub_ts_data$log10_timepoint[sub_ts_data$log10_timepoint=='-Inf']<-0
    sub_ts_means$log10_timepoint<-log10(sub_ts_means$timepoint)
    sub_ts_means$log10_timepoint[sub_ts_means$log10_timepoint=='-Inf']<-0
  }

  the_plot<-plot_cluster_traj(object,sub_ts_data,sub_ts_means,yaxis_name=yaxis_name,num_col=length(group_names)*num_clusters_per_row)
  return(the_plot)
}

#' @title Retrieve ancestor information from several TimeSeriesObjects
#'
#' @description This function queries multiple TimeSeriesObjects. Give a list of
#' ancestors, it checks which of its children GOs are found as significant in the
#' queried TimeSeriesObjects. Note that it can only do this for a singular timepoint
#' It then concatenates all of this information in a single dataframe. Of the extracted
#' information is a directionality, where an arrow is given (up or down) based on if the
#' identified 'child' is found to be higher expressed in the experiment (up) or in the
#' control (down).
#'
#' @param TS_object_paths A list of path to TimeSeriesObjects of interest
#' @param target_ancestors A vector of ancestors to query
#' @param target_timepoint The timepoint of interest
#' @param gpro_ont The gprofiler ontology to query
#'
#' @returns A data frame containing all the identified children of the ancestors
get_comparative_ancestor_data<-function(TS_object_paths,target_ancestors,target_timepoint,gpro_ont='GO:BP'){
  main_df<-data.frame(NULL)
  for(name in names(TS_object_paths)){
    load(TS_object_paths[[name]])
    gpro_res<-gprofiler_cluster_analysis(TS_object,gpro_ont,save_path = NULL)
    GO_clusters<-gpro_res[['GO_df']]
    #retain significant GOs only
    GO_clusters<-GO_clusters[GO_clusters$p_value<0.05,]

    ancestor_ont<-strsplit(gpro_ont,':')[[1]][2]
    GOs_ancestors_clust<-find_relation_to_ancestors(target_ancestors,GO_clusters,ontology =ancestor_ont)
    GOs_ancestors_clust$name<-rep(name,nrow(GOs_ancestors_clust))

    ts_data<-calculate_cluster_traj_data(TS_object,scale_feat=FALSE)
    mean_ts_data<-calculate_mean_cluster_traj(ts_data)

    if(target_timepoint!='mean'){
      ts_data<-ts_data[ts_data$timepoint==target_timepoint,]
      mean_ts_data<-mean_ts_data[mean_ts_data$timepoint==target_timepoint,]
    }
    GOs_ancestors_clust$direction<-rep('None',nrow(GOs_ancestors_clust))
    for(clust in unique(GOs_ancestors_clust$group_name)){
      temp_mean_df<-mean_ts_data[mean_ts_data$cluster==clust,]
      control_val<-temp_mean_df$trans_mean[temp_mean_df$group==TS_object@group_names[2]]
      exp_val<-temp_mean_df$trans_mean[temp_mean_df$group==TS_object@group_names[1]]
      if(target_timepoint=='mean'){
        control_val<-mean(control_val)
        exp_val<-mean(exp_val)
      }
      if(control_val>exp_val){
        direction<-'down'
      }else{
        direction<-'up'
      }
      num_of_timepoints<-length(unique(ts_data$timepoint))
      clust_size<-nrow(ts_data[ts_data$cluster==clust,])/2
      if(num_of_timepoints>1){
        clust_size<-clust_size/num_of_timepoints
      }
      new_clust_name<-paste0(clust,' (',clust_size,')')

      GOs_ancestors_clust$direction[GOs_ancestors_clust$group_name==clust]<-direction
      GOs_ancestors_clust$group_name[GOs_ancestors_clust$group_name==clust]<-new_clust_name
    }


    if(nrow(main_df)==0){
      main_df<-GOs_ancestors_clust
    }else{
      main_df<-rbind(main_df,GOs_ancestors_clust)
    }
  }

  main_df$name<-factor(main_df$name,levels=unique(main_df$name))

  return(main_df)
}

#' @title Plot the ancestor comparative dotplot
#'
#' @description Function which plots the data from the \code{get_comparative_ancestor_data}
#' function.
#'
#' @param plot_data The data as produced by /code{get_comparative_ancestor_data}
#' @param title The title to give to the plot. I reccommend including which
#' timepoint was used to ensure that there is no confusion in the interpretation.
#'
#' @import ggplot2 stringr
#'
#' @returns The ggplot2 dotplot
plot_comp_ancestor_plots<-function(plot_data,title){
  num_cols<-length(unique(plot_data$name))
  plt <-ggplot(plot_data, aes(x = group_name, y = term_name, color=ancestor_name,shape=direction)) +
    geom_point(aes(size=term_size,shape=direction),show.legend = TRUE)+
    scale_shape_manual(values=c("\u25B2","\u25BC"),breaks = c('up','down'))
  if(length(unique(plot_data$term_size))>1){
    plt<-plt+    scale_size(name   = "Term size",
                            breaks = c(min(plot_data$term_size),max(plot_data$term_size)),
                            labels = c(min(plot_data$term_size),max(plot_data$term_size)),
                            range = c(5,15))
  }
  plt<-plt+
    scale_y_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , "\n"),width = 60)) +
    theme_bw(base_size=15) +
    theme(axis.text.x=element_text(angle=-45,vjust=0)) +
    guides(shape = guide_legend(override.aes = list(size = 5)))+
    ylab("") +
    xlab("") +
    ggtitle(title)+
    facet_wrap(~name,scales='free_x',ncol=num_cols)

  return(plt)
}

