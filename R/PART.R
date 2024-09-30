# library(gprofiler2)
# library(tictoc)
# library('GOSemSim')

# packages_for_loading<-c('gprofiler2','tictoc','tibble','GOSemSim')
# suppressPackageStartupMessages(lapply(packages_for_loading, require, character.only = TRUE))

#' @title Format data for PART
#'
#' @description The function retrieves the counts from the time series object and filters
#' it to only contain relevant genes and samples. The data can also be scaled
#' if the parameter 'scale' is set to TRUE
#'
#'
#' @param object A time series object
#' @param target_genes Vector of genes to keep/use
#' @param scale Boolean indicating if the data should be scaled
#' @param target_samples Vector of samples to keep/use
#'
#' @return The timeseries object with the added PART count matrix
#'
#' @examples
#' TS_object<-create_example_object_for_R()
#' TS_object <- normalize_timeSeries_with_deseq2(time_object=TS_object)
#' #Perform conditional differential gene expression analysis
#' TS_object<-conditional_DE_wrapper(TS_object,vignette_run=TRUE)
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
#'
#' @export
#'
prep_counts_for_PART <-function(object,target_genes,scale,target_samples){
  #Check if the matrix already exists, if yes, end function
  if('part_matrix'%in% names(slot(object,'PART_results'))==TRUE){
    message('PART matrix already exists')
    return(object)
  }

  #Retrieve and prep assay
  cnts <- as.matrix(exp_matrix(object,'norm'))

  if (is.null(target_samples)==FALSE){
    cnts <- cnts[, target_samples]
  }

  target_genes <- target_genes[target_genes %in% rownames(cnts)]
  Y <- cnts[target_genes, ]
  if (scale) {
    Y <- t(scale(t(Y), center = TRUE, scale = TRUE))
  }

  object@PART_results[['part_matrix']] <- Y

  return (object)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @title compute PART clusters
#'
#' @description The function set's a seed for reproducibility of results, it then uses the
#' part function from the clusterGenomics package to establish which genes belong
#' to what clusters. The clusters are then ordered using hierarchical clustering.
#' This may result in some PART cluster being split, therefore the order is adjusted
#' to put them together. We first order using hierarchical clustering for visual purposes.
#'
#' @param object A timeseries object
#' @param part_recursion The number of recursions for PART calculation
#' @param part_min_clust The minimum number of genes per cluster
#' @param dist_param The distance parameter for clustering
#' @param hclust_param The hierarchical clustering method/parameter to be used
#' @param custom_seed The seed inputed (if any)
#' @param custom_matrix Allows the input of a custom matrix instead of taking it
#' from the object
#' @param return_as_object Boolean indicating if the results should be returned
#' within the submitted object or as a list
#' @param vignette_run Boolean indicating if this is for Vignettes, if so the function
#' will load the appropriate example data instead of performing the computation.
#'
#' @return The timeseries object with the PART results added
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
#'
#' @import tictoc
#' @importFrom tibble add_column
#' @importFrom stats hclust dist
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#'
#' @export
#'
compute_PART<-function(object,part_recursion=100,part_min_clust=10,
                       dist_param="euclidean", hclust_param="average",
                       custom_seed=NULL, custom_matrix=NULL,return_as_object=TRUE,
                       vignette_run=FALSE){
  if(vignette_run==TRUE){
    object<-part_load_results_vignettes(object)
    return(object)
  }

  #check if custom matrix was given
  if(is.null(custom_matrix)==TRUE){
    PART_res<-slot(object,'PART_results')
    #check if PART has already been calculated
    if('part_data' %in% names(PART_res)==TRUE){
      message('PART results already exists')
      return(object)
    }
    main_matrix<-PART_res$part_matrix
    message('computing PART clusters')

  }else{
    main_matrix<-custom_matrix
  }
  num_genes<-nrow(main_matrix)
  message(paste0('There are ',num_genes, ' in the matrix. If the computing time takes too long consider increasing the PART_l2fc_threshold.'))
  tic()#Start a timer for PART computation

  #Calculates the clustering using the 'part' algorithm
  calculated_clusters = part(main_matrix,B=part_recursion,minSize=part_min_clust,linkage=hclust_param)
  rowclust = hclust(dist(main_matrix,method=dist_param),method=hclust_param)

  clust_ordered <- unique(as.character(calculated_clusters$lab.hatK[rowclust$order]))

  #Rename clusters so they appear in order (i.e. 1,2,3,4)
  my_iter<-1
  clust_swap_df<-NULL
  for(val in unique(calculated_clusters$lab.hatK[rowclust$order])){
    if(is.null(clust_swap_df)==TRUE){
      clust_swap_df<-data.frame('original'=val,'replacement'=paste0('C',my_iter))
    }else{
      new_row<-c(val,paste0('C',my_iter))
      clust_swap_df<-rbind(clust_swap_df,new_row)
    }
    my_iter=my_iter+1
  }

  #Creates a save file for PART clusters
  part_data <- main_matrix[rowclust$order,]
  part_data <- as.data.frame(part_data)
  part_data <- add_column(part_data, gene_cluster = calculated_clusters$lab.hatK[rowclust$order], .before = 1)
  for(val in clust_swap_df$original){
    part_data$gene_cluster[part_data$gene_cluster==val]<-clust_swap_df$replacement[clust_swap_df$original==val]
  }
  #Re-order the data to fit the keep each cluster together
  part_data<-part_data[order(part_data$gene_cluster),]

  #Sets up the format for the color bar which will illustrate the different clusters found
  num_clusts <- length(unique(calculated_clusters$lab.hatK))
  cols <- brewer.pal(8, name = "Set3")[seq_len(min(8, num_clusts))]
  clust_cols <- colorRampPalette(colors = cols)(num_clusts)
  names(clust_cols) <- unique(part_data$gene_cluster)

  #Create the color vector to add to file
  col_vect<-clust_cols
  color_vector<-part_data$gene_cluster

  for (clust in unique(color_vector)){
    color_vector[color_vector==clust]<-col_vect[clust]
  }

  part_data <- add_column(part_data, cluster_col = color_vector, .before = 2)

  #Create cluster_map
  cluster_map <- part_data[,c(1,2)]
  colnames(cluster_map)=c('cluster','cluster_col')

  clust_info_list<-list(calculated_clusters=calculated_clusters,
                        clustered_rows=rowclust,
                        colored_clust_rows=clust_cols)

  PART_computation_time<-capture.output(toc())

  PART_params<-list(part_recursion=part_recursion, part_min_clust=part_min_clust,
                    dist_param=dist_param, hclust_param=hclust_param, custom_seed=custom_seed)

  if(return_as_object==TRUE){
    object@PART_results[['cluster_info']]<-clust_info_list
    object@PART_results[['part_data']]<-part_data
    object@PART_results[['cluster_map']]<-cluster_map
    object@PART_results[['comp_time']]<-PART_computation_time
    object@PART_results[['params_used']]<-PART_params

    return(object)
  }else{
    return_list<-list()
    return_list[['cluster_info']]<-clust_info_list
    return_list[['part_data']]<-part_data
    return_list[['cluster_map']]<-cluster_map
    return_list[['comp_time']]<-PART_computation_time
    return_list[['params_used']]<-PART_params
    return(return_list)
  }
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @title Gprofiler analysis
#'
#' @description Function which runs Gprofiler for PART clusters and stores the results
#' in the appropriate slot of the time object
#'
#' @param object A time series object
#' @param transcript_version_adjust str indicating if transcript version IDs have been submitted.
#' The string should reflect the pre-script of the IDs, for example, ENSG or ENSMUST.
#' If this is the case, the decimal and following number must be removed for the gprofiler analysis
#' as gprofiler recognizes genes, not transcripts.
#' @param gpro_sig whether all or only statistically significant results should be returned.
#' @param vignette_run Boolean indicating if this function is being run within
#' vignettes, if so it will bypass a network connection error and load mock data
#' otherwise the error will terminate the script and send the error back to the user.
#'
#' @return The updated object with the Gprofiler results
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
#' TS_object<-run_gprofiler_PART_clusters(TS_object,vignette_run=TRUE)
#'
#' @import gprofiler2
#' @import GOSemSim
#'
#' @export
run_gprofiler_PART_clusters<-function(object,transcript_version_adjust=NULL,gpro_sig=TRUE,vignette_run=FALSE){

  #Check if gprofiler results already exists
  if(length(slot(object,'Gprofiler_results'))>0){
    message('Gprofiler results already exist')
    return(object)
  }

  message('running Gprofiler on PART clusters')

  cmap<-slot(object,'PART_results')$cluster_map
  for (clust in unique(cmap$cluster)){
    gene_vect<-cmap[cmap==clust,]
    gene_vect<-as.vector(row.names(gene_vect))


    #Adjust for ENSG version
    if(is.null(transcript_version_adjust)==FALSE){
      ENSG_genes<-gene_vect[startsWith(gene_vect,transcript_version_adjust)]
      rem_decimal<-unlist(strsplit(ENSG_genes,'\\.'))[c(T,F)]
      #Remove transcripts from gene vect
      gene_vect<-gene_vect[!gene_vect %in% ENSG_genes]
      #Add the removed decimal version
      gene_vect<-c(gene_vect,rem_decimal)
    }

    #If function is called within Vignettes, mock data is used if the URL fails
    #This should only be done within the context of Vignettes
    if(vignette_run==TRUE){
      object<-part_gprofiler_vignettes(object)
      return(object)
    }else{#Vignette_run == FALSE
      message(paste0('Gprofiler for ',clust))
      gostres <- gost(query = gene_vect,organism = slot(object,'Gpro_org'),significant=gpro_sig)
    }
    if (is.null(gostres)==FALSE){
      object@Gprofiler_results[[clust]]<-gostres
    }
  }
  return(object)
}



#' @title Gprofiler analysis for Vignettes
#'
#' @description Function which loads gprofiler results for vignettes/examples
#'
#' @param object A time series object
#'
#' @return Updated object with gprofiler results
#' @export
part_gprofiler_vignettes<-function(object){
  object@Gprofiler_results<-PBMC_pre_loaded$gpro_res
  return(object)

}


#' @title PART analysis for Vignettes
#'
#' @description Function which loads PART results for Vignettes
#'
#' @param object A time series object
#'
#' @return Updated Object with PART results
#'
#'
#' @export
part_load_results_vignettes<-function(object){
  object@PART_results<-PBMC_pre_loaded$part_res
  return(object)
}
