# library(plotly)
# library(GO.db)
# library('GOSemSim')
# library(data.table)
# library(gprofiler2)
# library(htmlwidgets)
# library(dynamicTreeCut)
# library(stringr)

packages_for_loading<-c('plotly','GO.db','GOSemSim','data.table','gprofiler2',
                        'htmlwidgets','dynamicTreeCut','stringr')
suppressPackageStartupMessages(lapply(packages_for_loading, require, character.only = TRUE))

# pipeline/calculation functions  -------------------

#' @title Plot gprofiler results
#'
#' @description Creates a dataframe of GO's based on the gprofiler results stored in a time series
#' object. The dataframe is specific to the requested ontology
#'
#' If save_path is not null, the gprofiler results will be saved to the designated
#' location in both csv and html (interactive plot) format.
#'
#'
#' @param object A time series object
#' @param ontology the ontology that will be returned in the dataframe (ex: GO:BP)
#' @param save_path The folder path to save results if gprofiler results are to be saved
#' @param return_specific_cluster String to return the gost plot of a specific cluster
#'
#' @return GO_df A dataframe containing all the GOs (ID) found, their cluster, and the term name
#'
#' @export
#'
gprofiler_cluster_analysis<-function(object,ontology,save_path=NULL,return_specific_cluster=NULL){
  return_clust<-NULL
  #If a save_path was given, create folders to save the data
  if (is.null(save_path)==FALSE){
    dir.create(paste0(save_path,'/gprofiler_results'), showWarnings = FALSE)
    save_fig_location<-paste0(save_path,'/gprofiler_results/figures')
    save_data_location<-paste0(save_path,'/gprofiler_results/data_files')
    dir.create(save_fig_location, showWarnings = FALSE)
    dir.create(save_data_location, showWarnings = FALSE)
  }
  cluster_map<-object@PART_results$cluster_map
  GO_df<-as.data.frame(NULL)
  for (clust in names(object@Gprofiler_results)){
    gostres <- object@Gprofiler_results[[clust]]
    #store figure results
    p <- gostplot(gostres, capped = T, interactive = T)

    #Produce data results - conversion of parent column (collapse list)
    gprof_dta<-gostres$result
    for (idx in 1:nrow(gprof_dta)){
      gprof_dta$parents[idx]<-paste(unlist(gprof_dta$parents[idx]), collapse = '/')
    }
    gprof_dta$parents<-unlist(gprof_dta$parents)

    #Add to list for subsequent plotting
    if (ontology %in% gprof_dta$source == TRUE){
      extracted_results<-gprof_dta[c('term_id','term_name','p_value','term_size')][gprof_dta$source==ontology,]
      extracted_results$group_color<-rep(unique(cluster_map$cluster_col[cluster_map$cluster==clust]),nrow(extracted_results))
      extracted_results$group_name<-rep(unique(cluster_map$cluster[cluster_map$cluster==clust]),nrow(extracted_results))
      GO_df<-rbind(GO_df,extracted_results)
    }

    #If a save_path was given, save the data
    if (is.null(save_path)==FALSE){
      cluster_name<-strsplit(clust,'.csv')[[1]]
      save_name_fig<-paste0(save_fig_location,'/',cluster_name,'_overview.html')
      save_name_data<-paste0(save_data_location,'/',cluster_name,'_data.csv')
      #save plot
      saveWidget(p,file=save_name_fig)
      #If CORUM results, there may be some term names with commas
      gprof_dta$term_name<-gsub(",",".",gprof_dta$term_name)
      write.csv(x=gprof_dta,file=save_name_data,quote=FALSE)
      #Delete extra folder created by htmlwidgets
      unlink(paste0(save_fig_location,'/',cluster_name,'_overview_files'),recursive = T)

      if(is.null(return_specific_cluster)==F){
        if(cluster_name==return_specific_cluster){
          return_clust<-p
        }
      }
    }
  }
  GO_df<-as.data.frame(GO_df)
  GO_df[['-log10(padj)']]<--log10(p.adjust(GO_df$p_value, method = p.adjust.methods, n = length(GO_df$p_value)))
  return_list<-list(GO_df=GO_df,gost_clust=return_clust)
  return (return_list)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @title Calculate MDS matrix
#'
#' @description Function which calculates semantic distance and MDS from a dataframe of GOs
#'
#' The function uses a GOSemSim object to calculates the termdistance, followed
#' by a MDS calculation.
#' The resulting dataframe is formatted to be used in \code{plot_MDS}
#'
#' This function is used to create a GO term plot for WGCNA modules
#'
#'
#' @param GO_df The dataframe of GOs as returned by \code{gprofiler_cluster_analysis}
#' @param sem_data semantic similarity data as created by the godata function
#' @param measure the measuring method for semantic distance
#'
#' @return res.mds A dataframe containing the MDS dimensions, along with groupings and
#' colors which are equivalent to the cluster colors
#'
#' @export
#'
calculate_and_format_MDS<-function(GO_df,semantic_data,measure='Wang'){
  #Calculate term distance using Wang method
  TermDist=as.dist(1-mgoSim(GO_df$term_id,GO_df$term_id,semData=semantic_data,measure=measure,combine=NULL))
  # run MDS
  res.mds <-cmdscale(TermDist, eig =T, k = 2)

  # extract point values
  res.mds<-res.mds$points
  res.mds<-data.table(
    GO.ID=attr(res.mds,"dimnames")[[1]],
    Dim.1=res.mds[,1],
    Dim.2=res.mds[,2]
  )

  #Add module and term definitions
  res.mds$group_color <- factor(GO_df$group_color, levels=unique(GO_df$group_color) ,labels = unique(GO_df$group_color))
  res.mds$group_name<-GO_df$group_name
  res.mds$term<-GO_df$term_name


  return(res.mds)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @title merge_MDS modules
#'
#' @description Function finds and merges GOs from separate modules
#'
#' If a go appears as a part of two or more modules, the function will remove the duplicate
#' GOs and keep a single version which has all the relevant modules as it's source
#' ex: blue/yellow/brown if the GO was found in the blue, yellow, and brown cluster
#' It also finds the color which best mixes the colors of the cluster
#'
#'
#' @param the_data dataframe as returned by \code{calculate_and_format_MDS}
#'
#' @return the_data The updated inputted dataframe
#'
#' @export
#'
merge_duplicate_modules<-function(the_data){
  duplicated_GOs<-unique(the_data$GO.ID[duplicated(the_data$GO.ID)])

  dups<-the_data[the_data$GO.ID %in% duplicated_GOs,]
  list_new_groups<-c()
  my_new_df<-as.data.frame(NULL)

  for (dup_GO in unique(dups$GO.ID)){
    found_mods<-dups$group_name[dups$GO.ID==dup_GO]
    found_colors<-dups$group_color[dups$GO.ID==dup_GO]
    new_group<-paste(found_colors, collapse = "/")
    new_row<-dups[dups$GO.ID==dup_GO,]
    new_row$group_name<-paste(found_mods, collapse = "/")
    new_row<-new_row[!duplicated(new_row$GO.ID),]
    if (new_group!=''){
      new_row$group_color<-find_merged_color(new_group)
    }
    my_new_df<-rbind(my_new_df,new_row)
  }


  the_data<-the_data[!the_data$GO.ID %in% duplicated_GOs,]

  the_data<-rbind(the_data,my_new_df)

  return (the_data)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @title associate color to group for plotting
#'
#' @description Function which determines which color best represents color groupings (ex: blue/yellow)
#'
#' The function takes pairs of two colors, finds the 6th in between color using
#' coloRampPalette, and continues doing so until all colors in the group have
#' been included
#'
#'
#' @param group A vector of groups in the data
#'
#' @return new_col A vector of colors in the same order as the inputted group vector
#'
#' @examples
#'
#' @export
find_merged_color<-function(group){
  group_vect<-strsplit(group,'/')[[1]]
  idx=0
  new_col=NULL
  while (idx+1 != length(group_vect)){
    idx=idx+1
    if (is.null(new_col)==TRUE){
      pal<- colorRampPalette(c(group_vect[idx], group_vect[idx+1]))
      new_col<-pal(11)[6]
    }else{
      pal<- colorRampPalette(c(new_col, group_vect[idx+1]))
      new_col<-pal(11)[6]
    }
  }
  return(new_col)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @title Cluster GO terms using Wang measurement
#'
#' @description Function which clusters GO terms and returns the result in dataframe format
#'
#' This function uses the Wang measurement method with the ward.D2 aggregation method
#' clustering parameters are pre-set within the function.
#'
#' @param GO_df The dataframe of GOs as returned by \code{gprofiler_cluster_analysis}
#' @param semantic_data semantic similarity data as created by the godata function
#'
#' @return clust_df dataframe showing each GO with it's cluster
#'
#' @examples
#'
#' @export
find_clusters_from_termdist<-function(GO_df,semantic_data){

  dist=as.dist(1-mgoSim(GO_df$term_id,GO_df$term_id,semData=semantic_data,measure='Wang',combine=NULL))

  #Set-up pre-defined parameters for dendrogram
  Tree=list(
    tree=list(distance="Wang",aggreg.method="ward.D2",rotate=NULL),
    cut=list(static=NULL,dynamic=list(pamStage=TRUE,pamRespectsDendro=TRUE,deepSplit=2,minClusterSize =2)))

  ## aggregation
  dendro<-hclust(dist, method ='ward.D2')
  # create ordering vectored dendrogram
  ord<-order.dendrogram(as.dendrogram(dendro))

  # cut dynamic
  gp<-cutreeDynamic(dendro =dendro,method="hybrid",verbose=FALSE,
                    distM =as.matrix(dist),pamStage=Tree$cut$dynamic$pamStage,
                    pamRespectsDendro=Tree$cut$dynamic$pamRespectsDendro,
                    deepSplit=Tree$cut$dynamic$deepSplit,
                    minClusterSize =Tree$cut$dynamic$minClusterSize
  )
  # groups
  clust=unique(gp[ord])
  # new clusters names table
  clust=data.table(ini=clust,new=seq_along(clust))
  # convert gp to data.table
  gp<-data.table(ini=gp)
  # merge gp with new clusters names
  gp<-merge(gp,clust,by="ini",all.x=TRUE,sort=FALSE)
  # extract only new clusters names in the initial order
  gp<-gp$new
  # add go names
  names(gp)<-attr(dist,"Labels")

  #Converts the found clusters in a usable dataframe
  clust_df<-data.frame(names(gp)[order(unname(gp))],unname(gp)[order(unname(gp))])
  colnames(clust_df)=c('GO.ID','GO.cluster')

  for_merge<-GO_df[,c('term_id','group_color','group_name')]
  colnames(for_merge)<-c('GO.ID','GO.color','GO.clust.name')
  clust_df<-merge(clust_df,for_merge,by='GO.ID')
  clust_df<-clust_df[!duplicated(clust_df$GO.ID),]
  return(clust_df)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @title Find semantic similarity between clustered GOs
#'
#' @description Function which clusters GO terms and returns the result in dataframe format
#'
#' This function uses the Wang measurement method with the ward.D2 aggregation method
#' clustering parameters are pre-set within the function.
#'
#' Function was adapted from the calculate_SS function from the visEAGO package
#'
#' @param semantic_data semantic similarity data as created by the godata function
#' @param the_clusters GO terms with associated clutsers as created by \code{find_clusters_from_termdist}
#' @param selected_ont The ontology to be used
#' @param distance The distance algorithm to use (ex:'BMA')
#' @param measure The measuring algorithm to use (ex: 'Wang')
#'
#'
#' @return values A matrix of semantic distance for the clustered GOs
#'
#' @importFrom data.table ":="
#'
#' @examples
#'
#' @export
SS_GO_clusters<-function(semantic_data,the_clusters,selected_ont,distance,measure) {

  clusters<-the_clusters[,c('GO.ID','GO.cluster')]
  clusters<-as.data.table(clusters)
  # convert clusters in factor
  clusters[,"GO.cluster":=factor(clusters$GO.cluster,levels = unique(clusters$GO.cluster))]
  # get ancestors
  onto=switch(
    selected_ont,
    MF=GOMFANCESTOR,
    BP=GOBPANCESTOR,
    CC=GOCCANCESTOR
  )
  # convert in list
  onto=AnnotationDbi::as.list(onto)
  # keep enrich terms
  onto=stack(onto[clusters$GO.ID])

  # add enrich terms
  onto=rbind(onto,cbind(values=clusters$GO.ID,ind=clusters$GO.ID))

  # merge onto with GOclusters
  colnames(onto)<-c('values','GO.ID')
  onto<-merge(onto,clusters,by='GO.ID')
  onto<-onto[,c('GO.ID','GO.cluster','values')]
  onto<-onto[order(onto$GO.cluster,onto$GO.ID),]
  onto<-data.table(onto)

  # remove GO.ID
  onto=onto[,"GO.ID":=NULL]
  # count the term occurence by cluster
  onto=onto[,.N,by=c("GO.cluster","values")]
  # select the maxima term occurence by cluster
  onto=onto[onto[, .I[N==max(N)], by="GO.cluster"]$V1]
  # add IC to onto
  onto=data.table(onto,IC=semantic_data@IC[onto$values])
  # select the first maxima term IC by cluster
  onto=onto[onto[, .I[which.max(IC)], by="GO.cluster"]$V1]

  # extract GO terms
  GO<-data.table(
    AnnotationDbi::select(
      GO.db,
      keys=onto$values,
      column="TERM"
    )
  )

  # add sizes column to GOclusters
  onto=data.table(onto,GO)
  # convert clusters in factor
  clusters<-unstack(data.frame(clusters))
  # rename clusters
  names(clusters)<-paste(onto$GO.cluster,onto$GOID,onto$TERM,sep="_")

  ## compute clusters similarities

  # number of clusters
  n=length(clusters)

  # create default scores matrix
  scores <- matrix(NA, nrow = n, ncol = n)

  # create default scores matrix
  for (i in seq_along(clusters)){

    # extract the first cluster GO terms
    GO1 <-clusters[[i]]

    # extract second cluster GO terms
    for (j in seq_len(i)){

      # extract the second cluster GO terms
      GO2<-clusters[[j]]

      # calculate
      scores[i, j] <-mgoSim(
        GO1,
        GO2,
        semData = semantic_data,
        measure =measure,
        combine =distance
      )

      # add to upper part
      if (i != j) scores[j, i] <- scores[i, j]
    }
  }

  # add clusters name
  colnames(scores)<-names(clusters)
  row.names(scores)<-names(clusters)

  # return matrix
  values=list(stats::as.dist(1-scores))

  # add distance name to values
  names(values)=distance
  values<-as.matrix(values[[distance]])
  return(values)

}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @title Format clustered GO data for plotting
#'
#' @description Function which counts the number of GOs in each cluster and determines how many
#' GOs of each module appears in each cluster
#'
#' @param cluster_df A dataframe containing the GO terms, their cluster, and module(s)
#' as created by \code{find_clusters_from_termdist}
#'
#' @return cluster_nb_df A dataframe containing GO terms for the clusters, their associated
#' cluster, the number of GO terms in each cluster as well as the module dispertion
#' of the GOs within the cluster
#'
#' @export
#'
create_clustered_module_dataframe<-function(cluster_df){
  cluster_nb_df<-data.frame(matrix(ncol = 2, nrow = 0))
  for (clust in unique(cluster_df$GO.cluster[order(cluster_df$GO.cluster)])){
    new_row<-c(clust,length(cluster_df$GO.ID[cluster_df==clust]))
    cluster_nb_df<-rbind(cluster_nb_df,new_row)
  }
  colnames(cluster_nb_df)=c('GO.cluster','nb')
  cluster_nb_df$GO.cluster<-as.character(cluster_nb_df$GO.cluster)


  cluster_module_df<-cluster_df[order(cluster_df$GO.cluster),]
  merged_clusters<-c()
  for (clust in unique(cluster_module_df$GO.cluster)){
    cluster_dist<-table(cluster_module_df$GO.clust.name[cluster_module_df$GO.cluster==clust])
    #This is the ordering method for the appearance in plotly
    cluster_dist<-cluster_dist[order(-cluster_dist)]
    clust_string<-c()
    for (res in names(cluster_dist)){
      clust_string<-c(clust_string,paste0(res,': ',cluster_dist[[res]]))
    }
    clust_string<-paste(clust_string,collapse=", ")
    merged_clusters<-c(merged_clusters,clust_string)
  }
  #They are in the same order
  cluster_nb_df$GO.clust.name<-merged_clusters

  #Prepare the color_vector
  col_vect<-c()
  for(clust in cluster_nb_df$GO.clust.name){
    dominant_cluster<-strsplit(clust,split = ':')[[1]][1]
    col_vect<-c(col_vect,unique(cluster_df$GO.color[cluster_df$GO.clust.name==dominant_cluster]))
  }
  cluster_nb_df$GO.colors<-col_vect
  return(cluster_nb_df)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @title Ancestor querying for GOs
#'
#' @description Function which finds if a term is associated with a target ancestor
#'
#' The function will retrieve ancestry information from the annotationDbi of the
#' loaded library as well as the inputted ontology (BP by default). It will then
#' isolate and preserve the GO terms found to be related to one or several of the target
#' ancestors. If the GO is related to several ancestors, only the one highest in the
#' tree is linked to that GO.
#'
#' The function also associates a color to each ancestor, due to this, there is a
#' maximum of 9 ancestors that can be selected
#'
#' @param target_ancestors A vector of GO IDs to use as ancestors
#' @param GOs_to_check A dataframe created by either \code{gprofiler_for_conditional_DE}
#' or \code{gprofiler_cluster_analysis}
#' @param ontology The ontology to be used to search for ancestry
#'
#' @returns A updated dataframe where all non-ancestor related GOs were removed
#' GOs with a found ancestor are associated to it via a new column called 'ancestor'
#'
#' @export
find_relation_to_ancestors<-function(target_ancestors,GOs_to_check,ontology='BP'){
  # Get full list of GOs and ancestors
  # get ancestors
  onto=switch(
    ontology,
    MF=GOMFANCESTOR,
    BP=GOBPANCESTOR,
    CC=GOCCANCESTOR
  )
  onto=AnnotationDbi::as.list(onto)
  # Filter for GOs in our dataset
  onto=stack(onto[GOs_to_check$term_id])


  relations_found<-data.frame(term_id=NULL,ancestor=NULL)
  names_of_ancestor<-c()
  for (GO in unique(onto$ind)){
    sub_onto<-onto[onto$ind==GO,]
    if (any(sub_onto$values %in% target_ancestors)){
      found_ancestor<-sub_onto$values[sub_onto$values %in% target_ancestors][1]
      temp_df<-data.frame(term_id=GO,ancestor=found_ancestor)
      relations_found<-rbind(relations_found,temp_df)
    }
  }

  if(nrow(relations_found)==0){
    message('NO terms associated with given ancestors were found')
    return(data.frame(NULL))
  }

  #Give colors to each ancestor
  num_ancestors<-length(unique(relations_found$ancestor))
  cols <- RColorBrewer::brewer.pal(12, name = "Paired")
  names(cols)<-unique(relations_found$ancestor)
  relations_found$group_color<- cols[relations_found$ancestor]



  GOs_ancestor_found<-GOs_to_check[GOs_to_check$term_id %in% relations_found$term_id,]
  GOs_ancestor_found<-subset(GOs_ancestor_found, select=-c(group_color))
  GOs_ancestor_found<-merge(GOs_ancestor_found,relations_found,by='term_id')

  GOs_ancestor_found<-GOs_ancestor_found[order(GOs_ancestor_found$ancestor),]#Order by ancestor to ensure the same factorisation
  GOs_ancestor_found$group_color<-factor(GOs_ancestor_found$group_color,levels=unique(GOs_ancestor_found$group_color))

  # extract a named vector of all terms
  goterms <- Term(target_ancestors)
  goterms<-as.data.frame(goterms)
  goterms$ancestor<-row.names(goterms)
  colnames(goterms)<-c('ancestor_name','ancestor')
  GOs_ancestor_found<-merge(GOs_ancestor_found,goterms,by='ancestor')


  return(GOs_ancestor_found)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @title Read gprofiler data from results folder
#'
#' @description Function which reads the grprofiler data results from the time series analysis
#'
#' The function reads the information for each cluster and filters for the required
#' columns, ontology, and number of GOs. It then adds a -log10 FDR
#'
#' @param res_location The location of the gprofiler_results folder
#' @param ont The ontology to be filtered for
#' @param top_n The number of top GOs to be filtered for
#'
#' @return A dataframe with the requested information
#'
#' @export
read_gprofiler_results<-function(res_location='TS_results/',ont='REAC',top_n=NULL){
  GO_cluster<-data.frame(NULL)
  for(file in list.files(paste0(res_location,'gprofiler_results/data_files'))){
    read_path<-paste0(res_location,'gprofiler_results/data_files/',file)
    cluster<-strsplit(file,'_')[[1]][1]
    sub_df<-read.csv(read_path)
    sub_df<-sub_df[sub_df$source==ont,]
    sub_df$group_name<-rep(cluster,nrow(sub_df))

    sub_df<-sub_df[,c('term_id','term_name','p_value','source','term_size','group_name')]
    if(nrow(GO_cluster)==0){
      GO_cluster<-sub_df
    }else{
      GO_cluster<-rbind(GO_cluster,sub_df)
    }

  }
  if(nrow(GO_cluster)>0){
    if(is.null(top_n)==F){
      new_GO_df<-data.frame(NULL)
      for(clust in unique(GO_cluster$group_name)){
        sub_df<-GO_cluster[GO_cluster$group_name==clust,]
        if(nrow(sub_df)>top_n){
          sub_df<-sub_df[1:top_n,]
        }

        if(nrow(new_GO_df)==0){
          new_GO_df<-sub_df
        }else{
          new_GO_df<-rbind(new_GO_df,sub_df)
        }
      }
      GO_cluster<-new_GO_df
    }
    GO_cluster[['-log10(padj)']]<--log10(p.adjust(GO_cluster$p_value, method = p.adjust.methods, n = length(GO_cluster$p_value)))
    GO_cluster<-GO_cluster[,c('term_id','term_name','p_value','-log10(padj)','term_size','group_name')]
  }
  return(GO_cluster)
}

# plotting functions  -------------------

#' @title Plot nearest ancestor MDS (clustered GOs)
#'
#' @description Function which plots the MDS for clusters/nearest ancestor of the ancestor curation
#' dataframe
#'
#' @param the_data The dataframe as created by \code{calculate_and_format_MDS}
#'
#' @return p The plotly plot to enable either visualization or saving of the plot
#'
#' @export
#'
plot_ancestor_clust_MDS<-function(the_data,use_name=F){

  p<-plot_ly()
  color_ordered<-c()
  for (clust in sort(as.vector(unique(the_data$Ancestor)))){
    found_color<- as.vector(unique(the_data$group_color[the_data$Ancestor==clust]))
    color_ordered<-c(color_ordered,found_color)
  }
  # create trace
  p<-add_markers(
    p,
    data=the_data,
    x=~Dim.1,
    y=~Dim.2,
    color=~Ancestor,
    text = ~paste("Ancestor:",Ancestor,"<br>cluster(s):",group_name,"<br>GO.ID:",GO.ID,"<br>GO.name:",term,"<br>-log10(padj):",`-log10(padj)`),
    showlegend=TRUE,
    mode='markers',
    colors=color_ordered,
    sizes=c(20,50),
    size=~`-log10(padj)`,
    fill = ~'',#Doing this prevents a warning relating to line.width...somehow
    marker =list(
      sizemode = 'diameter',
      # color=~colors,
      # size =20,
      opacity = 0.4,
      showlegend=T
    )
  )
  # p
  return (p)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @title plot clustered MDS
#'
#' @description This function was adapted from the MDSplot() function from the visEAGO package
#'
#' @param main_matrix Matrix containing the data to be plotted as created by  \code{SS_GO_clusters}
#' @param cluster_info A dataframe with cluster names and the number of GOs in each cluster
#'
#' @return the plotly object for the MDS plot
#'
#' @importFrom data.table ":="
#'
#' @export
plot_clustered_mds<-function(main_matrix,cluster_info){
  # run MDS
  res.mds <-cmdscale(as.dist(main_matrix), eig = T, k = 2)

  # extract point values
  res.mds<-res.mds$points

  # convert to data.table
  res.mds<-data.table(
    GO.cluster=attr(res.mds,"dimnames")[[1]],
    Dim.1=res.mds[,1],
    Dim.2=res.mds[,2]
  )

  # custom text
  res.mds[,`:=`("text"=res.mds$GO.cluster,GO.cluster=gsub("_.+$","",res.mds$GO.cluster))]

  # add levels to GO.cluster
  res.mds$GO.cluster<-factor(
    res.mds$GO.cluster,
    levels=unique(res.mds$GO.cluster)
  )

  # add GO.ID for GO
  res.mds[,"text":=paste("GO.cluster:",text)]

  # add GO.ID for GO
  res.mds[,"text":=gsub("_GO","<br>GO.ID: GO",text)]

  # add GO.name for term definition
  res.mds[,"text":=gsub("_"," <br>GO.name: ",text)]

  ## plot

  # mds coordinates
  values=unlist(res.mds$Dim.1,res.mds$Dim.2)

  # init graph
  p<-plot_ly()
  # add count to res.mds
  res.mds=merge(
    res.mds,
    cluster_info,
    by="GO.cluster",
    sort=FALSE
  )
  # add count to text
  res.mds[,"text":=paste(sub("<br>.+$","",text),"<br>GO.count:",res.mds$nb,' – ',
                         cluster_info$GO.clust.name[as.numeric(res.mds$GO.cluster)],sub("^.+<br>GO.ID","<br>GO.ID",text))]

  # visualization uses color of dominant module
  colors<-c()
  for(clust in res.mds$GO.cluster){
    mod_string<-cluster_info$GO.colors[cluster_info$GO.cluster==clust]
    #Keep only first element (one with most GOs or alphabetical if match)
    colors<-c(colors,gsub(x=mod_string,pattern = ':.*',''))
  }

  # create trace
  p<-add_markers(
    p,
    data=res.mds,
    x=~Dim.1,
    y=~Dim.2,
    text=~text,
    showlegend=FALSE,
    sizes=c(20,50),
    size=~nb,
    fill = ~'',#Doing this prevents a warning relating to line.width...somehow
    marker =list(
      sizemode = 'diameter',
      opacity = 0.4,
      line=list(color=colors),
      color=colors
    ),
    visible=TRUE
  )

  # return the plot
  return(p)

}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @title Create dotplot for ancestor query results
#'
#' @description Function which creates a dotplot of all the GO terms found to be related
#' to selected ancestors. The function expects the dataframe from \code{find_relation_to_ancestors}
#'
#' The plot shows enirchment (up or down), to which ancestor a GO has been associated
#' as well as the size of the term. Each term is shown on the axis of the cluster
#' for wich it was found.
#'
#' @param ancestor_df The dataframe obtained from \code{find_relation_to_ancestors}
#' @param use_names Boolean indicating if the ancestor names or IDs should be used
#' @param enrichment_dta Boolean indicating if the plot has enrichment data
#'
#' @return The ggplot2 object for the created dotplot
#'
#' @export
dotplot_ancestors<-function(ancestor_df,use_names=F,enrichment_dta=T){

  if(use_names==T){
    ancestor_df$ancestor<-ancestor_df$ancestor_name
  }

  if (enrichment_dta==T){
    #Create shape vector and order based on enrichment
    #A certain shape is always associated to the same enrichment
    shape_vector<-c(22,24,25)
    names(shape_vector)<-c('both','up','down')
    shape_vector<-shape_vector[levels(ancestor_df$enrichment)]

    plt <-ggplot(ancestor_df, aes(x = timepoint, y = term_name, color=ancestor, fill= ancestor, shape=enrichment)) +
      geom_point(aes(color=ancestor,shape=enrichment,fill=ancestor)) +
      scale_shape_manual(values=unname(shape_vector)) +
      scale_fill_manual(values=levels(ancestor_df$group_color))+
      scale_colour_manual(values=levels(ancestor_df$group_color))+
      geom_point(aes(size=term_size),show.legend = T) +
      scale_size(name   = "GO size",
                 breaks = c(min(ancestor_df$term_size),max(ancestor_df$term_size)),
                 labels = c(min(ancestor_df$term_size),max(ancestor_df$term_size))) +
      scale_y_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , "\n"),width = 60))+
      theme_bw(base_size=15) +
      theme(axis.text.x=element_text(angle=-45,vjust=0)) +
      ylab("") +
      xlab("") +
      ggtitle("Ancestor analysis")
  }else{
    plt <-ggplot(ancestor_df, aes(x = group_name, y = term_name, color=ancestor, fill= ancestor)) +
      # geom_point(aes(color=ancestor,shape=ancestor,fill=ancestor)) +
      # scale_fill_manual(values=c("blue", "cyan4")) +
      scale_colour_manual(values=levels(ancestor_df$group_color))+
      geom_point(aes(size=`-log10(padj)`),show.legend = T) +
      # scale_size(name   = "GO size",
      #            breaks = c(min(ancestor_df$term_size),max(ancestor_df$term_size)),
      #            labels = c(min(ancestor_df$term_size),max(ancestor_df$term_size))) +
      scale_y_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , "\n"),width = 60))+
      theme_bw(base_size=15) +
      theme(axis.text.x=element_text(angle=-45,vjust=0)) +
      ylab("") +
      xlab("") +
      ggtitle("Ancestor analysis")
  }
  return(plt)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @title Create dotplot for gprofiler results
#'
#' @description Function which creates a dotplot for the provided gprofiler results
#'
#' The function creates the ggplot2 object for a dotplot which illustrates the GOs
#' of the give clusters.
#'
#' @param GO_clusters A dataframe with the group_name(cluster), term_name, -log10FDR,
#' and the term_size
#' @param ont The ontology used, it will be added to the title
#' @param top_n The number of top GOs used for the filter, will be added to the title
#'
#' @return Returns the ggplot2 object for the plot
#'
#' @export
custom_gpro_dotplot<-function(GO_clusters,ont,top_n){
  plt <-ggplot(GO_clusters, aes(x = group_name, y = term_name, color=`-log10(padj)`)) +
    scale_color_gradient(low = "blue", high = "red") +
    geom_point(aes(size=term_size),show.legend = T)
  if(length(unique(GO_clusters$term_size))>1){
    plt<-plt+    scale_size(name   = "Term size",
                            breaks = c(min(GO_clusters$term_size),max(GO_clusters$term_size)),
                            labels = c(min(GO_clusters$term_size),max(GO_clusters$term_size)))
  }
  plt<-plt+
    scale_y_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , "\n"),width = 60)) +
    theme_bw(base_size=15) +
    theme(axis.text.x=element_text(angle=-45,vjust=0)) +
    ylab("") +
    xlab("") +
    ggtitle(paste0(ont,", top ",top_n))

  return(plt)
}

#' @title Create MDS plot
#'
#' @description The function will plot the GOs for each module while preserving their module name
#' as their color.
#'
#' @param the_data The dataframe as created by \code{calculate_and_format_MDS}
#'
#' @return p The plotly plot to enable either visualization or saving of the plot
#'
#'
#' @export
plot_MDS<-function(the_data){
  p<-plot_ly()
  color_ordered<-c()
  for (clust in sort(as.vector(unique(the_data$group_name)))){
    found_color<- as.vector(unique(the_data$group_color[the_data$group_name==clust]))
    color_ordered<-c(color_ordered,found_color)
  }
  # create trace
  p<-add_markers(
    p,
    data=the_data,
    x=~Dim.1,
    y=~Dim.2,
    color=~group_name,
    text = ~paste("modules/clusters:",group_name,"<br>GO.ID:",GO.ID,"<br>GO.name:",term),
    showlegend=TRUE,
    mode='markers',
    colors=color_ordered,
    marker =list(
      # color=~colors,
      size =20,
      opacity = 0.4,
      showlegend=T
    )
  )
  # p
  return (p)
}

# Plotting wrapper functions  -------------------

#' @title Create gprofiler dotplots
#'
#' @description Wrapper function for the creation of a dotplot which summarises the gprofiler findings
#' for a specified ontology.
#'
#' The function reads the results for a location. It will automatically search for the
#' 'gprofiler_results/data_files' folder in the procided location. This is should be located
#' in the TS_results folder, but can change based on the user.
#'
#' The function will attempt to create the plot with proper width and height,
#' but these can be over-written by providing them to the function.
#'
#' The function will save the plot in the location provided through the
#' gpro_file_location parameter.
#'
#' @param gpro_file_location The location where the 'gprofiler_results' folder is
#' @param target_ontology The targeted ontology ex: 'REAC' or 'GO:BP'
#' @param top_n The number of top GOs to plot per cluster
#' @param custom_width A custom value for the width of the plot
#' @param custom_height A custom value for the height of the plot
#' @param return_plot Boolean indicating if the plot should be returned
#'
#' @return if specified, the function will return the ggplot2 object for the dotplot
#'
#' @export
GO_dotplot_wrapper<-function(gpro_file_location,target_ontology,top_n,custom_width=NULL,custom_height=NULL,return_plot=F){
  GO_top_cluster<-read_gprofiler_results(gpro_file_location,target_ontology,top_n)
  if(nrow(GO_top_cluster)>0){
    GO_top_cluster <- GO_top_cluster[order(GO_top_cluster[,'term_id'],-GO_top_cluster[,'-log10(padj)']),]
    # GO_top_cluster <- GO_top_cluster[!duplicated(GO_top_cluster$term_id),]
    plt<-custom_gpro_dotplot(GO_top_cluster,target_ontology,top_n)
    if(is.null(custom_height)==T){
      custom_height<-3+nrow(GO_top_cluster)*0.25
    }
    if(is.null(custom_width)==T){
      custom_width<-6+(length(unique(GO_top_cluster$group_name))*0.30)
    }
    svg(paste0(gpro_file_location,'/gprofiler_results/',target_ontology,'_dotplot.svg'),width=custom_width,height=custom_height)
    print(plt)
    trash<-capture.output(dev.off())

    if(return_plot==T){
      return(plt)
    }
  }
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @title Create MDS and nearest ancestor MDS plots
#'
#' @description Wrapper function which makes the necessary calls to plot a multi dimensional
#' scaling (MDS) plot for all of the terms found in each cluster (for a specified ontology).
#' The function calculates and plots a nearest ancestor version of the MDS using the
#' BMA and WANG approach. Each term is brought up to it's nearest ancestor.
#'
#' Each plot is saved as an interactive plot in TS_results
#'
#' @param GO_clusters The dataframe of GOs as returned by \code{gprofiler_cluster_analysis}
#' @param sem_data semantic similarity data as created by the godata function
#' @param sem_ontology The ontology that was used to calculate the semantic data
#' and that will be plotted
#' @param target_dir string indicating the location where the results will be saved
#' @param return_plot boolean indicating if the plots should be returned
#'
#' @return if specified, will return a list containing the two plotly objects for
#' both the GO term MDS and clustered GO term MDS.
#'
#' @export
wrapper_MDS_and_MDS_clusters<-function(GO_clusters,sem_data,sem_ontology,target_dir='TS_results/gprofiler_results/',return_plot=F){
  plot_data<-calculate_and_format_MDS(GO_clusters,sem_data)
  plot_data<-merge_duplicate_modules(plot_data)
  my_plot<-plot_MDS(plot_data)

  htmlwidgets::saveWidget(as_widget(my_plot), paste0(target_dir,"MDS_GO_terms.html"))
  unlink(paste0(target_dir,"MDS_GO_terms_files"),recursive = T)

  found_clusters<-find_clusters_from_termdist(GO_clusters,sem_data)
  # calculate semantic similarites between clusters of GO terms
  calculated_SS<-SS_GO_clusters(sem_data,found_clusters,sem_ontology,distance="BMA",measure='Wang')
  clustered_module_df<-create_clustered_module_dataframe(found_clusters)

  clust_plot<-plot_clustered_mds(calculated_SS,clustered_module_df)
  htmlwidgets::saveWidget(as_widget(clust_plot), paste0(target_dir,"MDS_GO_clusters.html"))
  unlink(paste0(target_dir,"MDS_GO_clusters_files"),recursive = T)

  if(return_plot==T){
    return(list(MDS_term=my_plot,MDS_clust=clust_plot))
  }
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @title Create ancestor query plots
#'
#' @description Wrapper function which makes the necessary calls to plots ancestor curation plots
#'
#' These plots find all terms which are children of the provided ancestors within the
#' specified ontology.
#'
#' The function will plot and save a dotplot as well as an interactive MDS plot
#'
#' @param GO_df The dataframe of GOs as returned by \code{find_relation_to_ancestors}
#' @param sem_data semantic similarity data as created by the godata function
#' @param use_names boolean indicating if names or IDs should be used
#' @param target_dir string indicating the save location of the plots
#' @param return_plot boolean indicating if the plot should be returned
#'
#' @return if specified, will return a list containing the ggplot2 object for the dotplot
#' and the plotly object for the MDS plot.
#'
#' @export
wrapper_ancestor_curation_plots<-function(GO_df,sem_data,use_names=T,target_dir='TS_results/',return_plot=F){

  if(is.null(GO_df)==T){
    return(NULL)
  }else if(nrow(GO_df)>0){
    dir.create(paste0(target_dir,'ancestor_plots'))

    custom_width<-8+(length(unique(GO_df$group_name))*0.30)
    custom_height<-2+(nrow(GO_df)*0.25)

    my_dotplot<-dotplot_ancestors(GO_df,enrichment_dta=F,use_names=use_names)

    svg(paste0(target_dir,'ancestor_plots/ancestors_clusters.svg'),width=custom_width,height=custom_height)
    print(my_dotplot)
    dev.off()
  }else{
    my_dotplot<-NULL
  }


  if(length(unique(GO_df$term_id))>1 & length(unique(GO_df$ancestor_name))>1){
    plot_data<-calculate_and_format_MDS(GO_df,sem_data)
    plot_data<-merge_duplicate_modules(plot_data)

    if(use_names==T){
      for_merger<-GO_df[,c('term_id','ancestor_name','-log10(padj)')]
    }else{
      for_merger<-GO_df[,c('term_id','ancestor','-log10(padj)')]
    }

    colnames(for_merger)=c('GO.ID','Ancestor','-log10(padj)')
    plot_data<-merge(plot_data,for_merger,by='GO.ID')

    my_MDS<-plot_ancestor_clust_MDS(plot_data)
    htmlwidgets::saveWidget(as_widget(my_MDS), paste0(target_dir,"ancestor_plots/cluster_ancestors.html"))
    unlink(paste0(target_dir,"ancestor_plots/cluster_ancestors_files"),recursive = T)

    write.csv(GO_df,paste0(target_dir,'ancestor_plots/ancestor_data.csv'),row.names=F)
  }else{
    my_MDS<-NULL
  }

  if(return_plot==T){
    return(list(dotplot=my_dotplot,MDS=my_MDS))
  }

}
