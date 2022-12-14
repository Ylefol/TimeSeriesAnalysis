# library(clusterGenomics)
# library(gprofiler2)
# library(tictoc)
# library('GOSemSim')

packages_for_loading<-c('clusterGenomics','gprofiler2','tictoc','tibble','GOSemSim')
suppressPackageStartupMessages(lapply(packages_for_loading, require, character.only = TRUE))

#' @title Format data for PART
#'
#' @description The function retrieves the counts from the time series object and filters
#' it to only contain relevant genes and samples. The data can also be scaled
#' if the parameter 'scale' is set to T
#'
#'
#' @param object A time series object
#' @param target_genes Vector of genes to keep/use
#' @param scale Boolean indicating if the data should be scaled
#' @param target_samples Vector of samples to keep/use
#'
#' @return The timeseries object with the added PART count matrix
#'
#' @export
#'
prep_counts_for_PART <-function(object,target_genes,scale,target_samples){
  #Check if the matrix already exists, if yes, end function
  if('part_matrix'%in% names(object@PART_results)==T){
    message('PART matrix already exists')
    return(object)
  }

  #Retrieve and prep assay
  cnts <- object@count_matrix$norm

  if (is.null(target_samples)==F){
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
#' @description The function set's a see for reproducibility of results, it then uses the
#' part function from the clusterGenomics package to establish which genes belong
#' to what clusters.
#'
#' @param object A timeseries object
#' @param save_name the file name to save the PART data, if NULL no saving
#' @param part_recursion The number of recursions for PART calculation
#' @param part_min_clust The minimum number of genes per cluster
#' @param dist_param The distance parameter for clustering
#' @param hclust_param The hierarchical clustering method/parameter to be used
#' @param custom_seed An integer to be used to set a custom seed
#' which gives the option of reproducibility from the PART calculation
#' @param custom_matrix Allows the input of a custom matrix instead of taking it
#' from the object
#' @param return_as_object Boolean indicating if the results should be returned
#' within the submitted object or as a list
#'
#' @return The timeseries object with the PART results added
#'
#' @export
#'
compute_PART<-function(object,part_recursion=100,part_min_clust=10,
                           dist_param="euclidean", hclust_param="average",
                           custom_seed=NULL,custom_matrix=NULL,return_as_object=T){

  #check if custom matrix was given
  if(is.null(custom_matrix)==T){
    #check if PART has already been calculated
    if('part_data' %in% names(object@PART_results)==T){
      message('PART results already exists')
      return(object)
    }
    main_matrix<-object@PART_results$part_matrix
    message('computing PART clusters')

  }else{
    main_matrix<-custom_matrix
  }

  #Sets a seed for reproducibility
  if (is.null(custom_seed)==F){
    set.seed(as.character(custom_seed))
  }

  tic()#Start a timer for PART computation

  #Calculates the clustering using the 'part' algorithm
  calculated_clusters = part(main_matrix,B=part_recursion,minSize=part_min_clust)
  rowclust = hclust(dist(main_matrix,method=dist_param),method=hclust_param)

  clust_ordered <- unique(as.character(calculated_clusters$lab.hatK[rowclust$order]))

  #Rename clusters so they appear in order (i.e. 1,2,3,4)
  my_iter<-1
  new_vect<-c()
  for(val in unique(calculated_clusters$lab.hatK[rowclust$order])){
    rep_value<-unname(table(calculated_clusters$lab.hatK[rowclust$order])[as.character(val)])
    new_vect<-c(new_vect,rep(as.character(my_iter),rep_value))
    my_iter<-my_iter+1
  }

  #Sets up the format for the color bar which will illustrate the different clusters found
  num_clusts <- length(unique(calculated_clusters$lab.hatK))
  cols <- RColorBrewer::brewer.pal(8, name = "Set3")[seq_len(min(8, num_clusts))]
  clust_cols <- grDevices::colorRampPalette(colors = cols)(num_clusts)
  names(clust_cols) <- unique(as.character(calculated_clusters$lab.hatK))

  #Create the color vector to add to file
  col_vect<-clust_cols[clust_ordered]
  names(col_vect)=seq(1,length(col_vect))
  color_vector<-c()

  #Creates a save file for PART clusters
  part_data <- main_matrix[rowclust$order,]
  part_data <- as.data.frame(part_data)
  part_data <- add_column(part_data, gene_cluster = new_vect, .before = 1)
  for (clust in unique(part_data$gene_cluster)){
    color_vector<-c(color_vector,rep(col_vect[clust],nrow(part_data[part_data$gene_cluster==clust,])))
  }
  part_data <- add_column(part_data, cluster_col = color_vector, .after = 1)

  #Create cluster_map
  cluster_map <- part_data[,c(1,2)]
  cluster_map$gene_cluster <- paste("C", cluster_map$gene_cluster,sep='')
  colnames(cluster_map)=c('cluster','cluster_col')

  clust_info_list<-list(calculated_clusters=calculated_clusters,
                        clustered_rows=rowclust,
                        colored_clust_rows=clust_cols)

  PART_computation_time<-capture.output(toc())

  PART_params<-list(part_recursion=part_recursion, part_min_clust=part_min_clust,
                    dist_param=dist_param, hclust_param=hclust_param,
                    custom_seed=custom_seed)

  if(return_as_object==T){
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
#'
#' @return The updated object with the Gprofiler results
#'
#' @export
run_gprofiler_PART_clusters<-function(object){

  #Check if gprofiler results already exists
  if(length(object@Gprofiler_results)>0){
    message('Gprofiler results already exist')
    return(object)
  }

  message('running Gprofiler on PART clusters')

  cmap<-object@PART_results$cluster_map
  for (clust in unique(cmap$cluster)){
    message(paste0('Gprofiler for ',clust))
    gene_vect<-cmap[cmap==clust,]
    gene_vect<-as.vector(row.names(gene_vect))
    gostres <- gost(query = gene_vect,organism = object@Gpro_org)
    if (is.null(gostres)==F){
      object@Gprofiler_results[[clust]]<-gostres
    }
  }
  return(object)
}
