

# library(ComplexHeatmap)
# library(ggplot2)
# library(ggrepel)
# library(grid)
# library(DESeq2, include.only = c("vst","plotPCA"))
# library(limma, include.only = c("plotMDS"))
# library(reshape2)

# packages_for_loading<-c('DESeq2','limma','grid','ggrepel','ggplot2','ComplexHeatmap','reshape2','stringr','stringi')
# suppressPackageStartupMessages(lapply(packages_for_loading, require, character.only = TRUE))

# wrapper for result creation -------------------

#' @title Creates differential expression results
#'
#' @description The function will create a folder called DE_results_(conditonal or temporal),
#' it will then create a summary heatmap of all the differential expression experiments
#' done in that type.
#'
#' The function will then iterate over the differential gene expression experiments done in
#' either DE_type (conditional or temporal) and creates two csv files, one containing
#' the raw differential expression data and the other is a subset of the raw data
#' containing only the significant genes
#' The function also creates a volcano plot and MA plot for each experiment
#'
#' The individual differential expression results are stored in folders which
#' carry the name of the experiment
#'
#' @param object A timeseries object
#' @param DE_type The type of differential gene expression done that the function
#' will use to create results. Either conditional or temporal
#' @param genes_of_interest A vector of genes of interest, can be empty (c())
#' @param results_folder The main folder where the results will be stored
#' @param adjust_missing_temp_samples Boolean indicating if missing temporal samples should be
#' compensated for (NA introduced). Currently, there is no reason to have this set to
#' False.
#' @param do_SVGs Boolean indicating if SVGs should be saved for the heatmaps
#'
#' @import ggplot2
#'
#' @return None
#'
#' @export
plot_wrapper_DE_results<-function(object,DE_type,genes_of_interest=c(),results_folder='TS_results/',adjust_missing_temp_samples=TRUE,do_SVGs=TRUE){
  #Create path to DE type and create directory
  main_path<-paste0(results_folder,'DE_results_',DE_type,'/')
  dir.create(main_path)

  message('Creating summary heat map')
  #Create a overview heatmap of the differential expression results for that type
  custom_heatmap_wrapper(object,DE_type=DE_type,log_transform=TRUE,
                         plot_file_name = paste0(main_path,'custom_heat_',DE_type),
                         adjust_missing_temp_samples,do_SVGs=do_SVGs)

  DE_res<-slot(object,'DE_results')[[DE_type]]

  #Iterate over the experiments in the DE_type and create csv results, volcano plot, and MA plots
  for (DE_name in names(DE_res)){
    message(paste0('creating results for ',DE_name))
    #Create a save path and the associated folder
    save_path<-paste0(main_path,DE_name,'/')
    dir.create(save_path)

    #Extract filter parameters for significance
    DE_p_filter<-slot(object,'DE_p_filter')
    DE_p_thresh<-slot(object,'DE_p_thresh')
    DE_l2fc_thresh<-slot(object,'DE_l2fc_thresh')

    #Create csv results for raw and significant differential expression data
    create_DE_data_results(object,DE_type=DE_type,exp_name=DE_name,
                           save_location=save_path)

    #Extract raw differential expression data
    DE_results<-DE_res[[DE_name]][['DE_raw_data']]

    #Create and save volcano plot
    v_plot<-volcanoplot_alt(DE_results, genes_of_interest=genes_of_interest,
                            filter_choice=DE_p_filter,l2FC_thresh=DE_l2fc_thresh,
                            p_thresh=DE_p_thresh, plot_title=DE_name, label_top_n=0,show_non_sig_interest=FALSE)
    ggsave(paste0(save_path,"volcano_plot.png"),dpi=300,width=21, height=19, units='cm',
           plot=v_plot)

    #Create and save the MA plot
    #No genes of interest, no labeling wanted
    ma_plot<-maplot_alt(DE_results, genes_of_interest=c(),
                        filter_choice=DE_p_filter,l2FC_thresh=DE_l2fc_thresh,
                        p_thresh=DE_p_thresh, plot_title=DE_name)
    ggsave(paste0(save_path,"MA_plot.png"),dpi=300,width=21, height=19, units='cm',
           plot=ma_plot)


    #REMOVED TEMPORARLY
    # #Plot and save the PCA plot
    # pca_plot<-plot_PCA_TS(time_object=object,exp_name=DE_name,DE_type=DE_type)
    # ggsave(paste0(save_path,"PCA_plot.png"),dpi=300,width=21, height=19, units='cm',
    #        plot=pca_plot)
  }
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @title Creates differential expression results for Vignettes
#'
#' @description The function will generate two summary heatmaps of the differential
#' gene expression results. The function was originally written to be used in Vignettes.
#' Full sets of differential gene expression results can be obtained with the
#' \code{plot_wrapper_DE_results} function
#'
#' @param object A timeseries object
#'
#' @import ggplot2
#'
#' @return list of heatmaps (contains conditional and temporal plot)
#'
#' @examples
#' TS_object<-create_example_object_for_R()
#' TS_object <- normalize_timeSeries_with_deseq2(time_object=TS_object)
#' #Perform conditional differential gene expression analysis
#' TS_object<-conditional_DE_wrapper(TS_object)
#' TS_object<-temporal_DE_wrapper(TS_object)
#' plot_list<-DE_plots_vignettes(TS_object)
#' @export
DE_plots_vignettes<-function(object){

  cond_heat_plot<- custom_heatmap_wrapper(object,DE_type='conditional',log_transform=TRUE,
                                          plot_file_name = NULL,do_SVGs=FALSE)

  temp_heat_plot<- custom_heatmap_wrapper(object,DE_type='temporal',log_transform=TRUE,
                                          plot_file_name = NULL,do_SVGs=FALSE)

  heat_list<-list('conditional'=cond_heat_plot,'temporal'=temp_heat_plot)

  return(heat_list)

}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @title Plot cluster trajectory for all PART clusters
#'
#' @description Function which plots the trajectory of all the clusters given
#' The function will split the clusters into groups of 8 to keep
#' the figures clean. If multiple figures are created, they will be names
#' 1_of_x until it reaches x_of_x.
#'
#' @param object A timeseries object
#' @param cluster_traj_dta The trajectory data for all clusters being calculated
#' The data is calculated/obtained from \code{calculate_cluster_traj_data} function
#' @param mean_cluster_traj_dta The mean value for each clusters trajectory,
#' the data is calculated/obtained from \code{calculated_mean_cluster_traj} function
#' @param log_TP Boolean indicating if timepoints (x axis) should be log10 transformed
#' @param plot_name The name given to the plot file as it is saved
#'
#' @return None
#'
#' @export
wrapper_cluster_trajectory<-function(object,cluster_traj_dta,mean_cluster_traj_dta,log_TP=FALSE,plot_name='Ctraj'){
  clust_order<-unique(cluster_traj_dta[,c('cluster','nGenes')])
  clust_order<-clust_order$cluster[order(-clust_order$nGenes)]
  num_needed_figures<-ceiling(length(clust_order)/8)
  #Iterate over number of necessary figures
  for (idx in 1:num_needed_figures){
    #Filter the cluster_map dataframe for the required clusters
    max_clust=8*idx
    if (idx==1){
      min_clust<-1
    }else{
      min_clust<-8*(idx-1)+1
    }
    clusters_to_plot<-clust_order[min_clust:max_clust]
    clusters_to_plot<-clusters_to_plot[!is.na(clusters_to_plot)]#Remove NAs

    sub_ts_data<-cluster_traj_dta[cluster_traj_dta$cluster %in% clusters_to_plot,]
    sub_ts_data<-sub_ts_data[order(match(sub_ts_data$cluster,clusters_to_plot)),]

    sub_ts_means<-mean_cluster_traj_dta[mean_cluster_traj_dta$cluster %in% clusters_to_plot,]
    sub_ts_means<-sub_ts_means[order(match(sub_ts_means$cluster,clusters_to_plot)),]

    if (num_needed_figures > 1){
      save_name<-paste0(plot_name,'_',idx,'_of_',num_needed_figures,'.svg')
    }else{
      save_name<-paste0(plot_name,'.svg')
    }

    sub_ts_data$labels<-factor(sub_ts_data$labels,levels = unique(sub_ts_data$labels))
    sub_ts_means$labels<-factor(sub_ts_means$labels,levels = unique(sub_ts_means$labels))

    cluster_num<-length(clusters_to_plot)
    number_rows<-ceiling(cluster_num/2)
    custom_height<-3*number_rows
    if (cluster_num==1){
      custom_width<-6
    }else{
      custom_width<-12
    }

    if(log_TP==TRUE){
      sub_ts_data$log10_timepoint<-log10(sub_ts_data$timepoint)
      sub_ts_data$log10_timepoint[sub_ts_data$log10_timepoint=='-Inf']<-0
      sub_ts_means$log10_timepoint<-log10(sub_ts_means$timepoint)
      sub_ts_means$log10_timepoint[sub_ts_means$log10_timepoint=='-Inf']<-0

    }

    the_plot<-plot_cluster_traj(object,sub_ts_data,sub_ts_means)
    svg(save_name,width=custom_width,height=custom_height)
    print(the_plot)
    dev.off()
  }
}

# Differential gene expression plots -------------------



#' @title Create a volcano plot
#'
#' @description The function creates a 'four quadrant' volcano plot, where the FDR and log2FoldChange
#' thresholds dictate the significant up-regulated category, the significant downregulated
#' category, the significant low regulation category and the non-significant category
#'
#' Genes of interest are labeled in a rectangle for visibility and will have the same
#' color as the category in which they are in.
#' Top significant genes are also labelled in rectangles but will have black text
#' in order to distinguish them from the genes of interest
#'
#' The plot is created using ggplot2, to save a plot the ggsave() function is recommended
#' It is also recommended to use the following parameters to save the plot.
#' dpi=300 width=21 height=19 units='cm'
#'
#' @param DE_res The differential expression results to be plotted
#' @param genes_of_interest A vector containing gene names to be labelled in the plot
#' To not label any genes, leave as default or provide an empty vector.
#' @param filter_choice Either padj or pvalue, the choice will be used to filter for
#' significance
#' @param l2FC_thresh The log2FoldChange threshold used to establish significance
#' @param p_thresh The pvalue or padj threshold used to establish significance
#' @param plot_title The title to be give to the plot
#' @param label_top_n The number of genes to label. Genes will be taken in order of
#' significance where the genes with the lowest adjusted p-value are taken first.
#' @param show_non_sig_interest A boolean to indicate if the non-significant genes of
#' interest should be shown.
#'
#' @return The ggplot object for the volcano plot
#'
#' @examples
#' TS_object<-create_example_object_for_R()
#' TS_object <- normalize_timeSeries_with_deseq2(time_object=TS_object)
#' #Perform conditional differential gene expression analysis
#' TS_object<-conditional_DE_wrapper(TS_object)
#' DE_res<-slot(TS_object,'DE_results')$conditional$IgM_vs_LPS_TP_1$DE_raw_data
#' v_plot<-volcanoplot_alt(DE_res = DE_res)
#'
#' @import ggplot2
#' @importFrom  ggrepel geom_label_repel
#'
#' @export
volcanoplot_alt <- function(DE_res,genes_of_interest=c(),filter_choice='padj',l2FC_thresh=1,
                            p_thresh=0.05,plot_title="Volcano plot",label_top_n=0,show_non_sig_interest=TRUE){

  DE_res$Significance <- NA
  sigUp <- which(DE_res$log2FoldChange>l2FC_thresh & DE_res[filter_choice]<p_thresh)
  DE_res[ sigUp, "Significance"] <- paste0("up-reg with ",filter_choice,"<",p_thresh," (",length(sigUp),")")

  #subset and fill columns based on significant downregulation
  sigDown <- which(DE_res$log2FoldChange<(-l2FC_thresh) & DE_res[filter_choice]<p_thresh)
  DE_res[ sigDown, "Significance"] <- paste0("down-reg with ",filter_choice,"<",p_thresh," (",length(sigDown),")")

  #subset and fill columns based on non-significant upregulation
  nonSig <- which(DE_res[filter_choice]>=p_thresh)
  DE_res[ nonSig, "Significance"] <- paste0("non-significant with ",filter_choice,">",p_thresh," (",length(nonSig),")")

  lowSig <- which(abs(DE_res$log2FoldChange)<=l2FC_thresh & DE_res[filter_choice]<p_thresh)
  DE_res[ lowSig, "Significance"] <- paste0("low-regulation with ",filter_choice,"<",p_thresh," (",length(lowSig),")")

  #Find the lowest significant point to draw a horizontal ax-line at that point on the axis
  #Sort using the filter choice
  DE_res <- DE_res[order(DE_res[[filter_choice]]),]
  #Find the 'first' non significant gene
  gene_for_sig_line<-DE_res$gene_id[DE_res$Significance==paste0("non-significant with ",filter_choice,">",p_thresh," (",length(nonSig),")")][1]
  #Create a dataframe of necessary columns and convert pvalue to -log10 (as is it's form in the plot)
  find_lowest_sig<-DE_res[,c('gene_id','pvalue')]
  find_lowest_sig$pvalue<- -log10(find_lowest_sig$pvalue)
  #Extract value based on the gene.
  lowest_sig<-find_lowest_sig$pvalue[find_lowest_sig$gene_id==gene_for_sig_line]

  labs_data <- DE_res[DE_res$gene_id %in% genes_of_interest, ]
  if (show_non_sig_interest==FALSE){
    labs_data <- labs_data[labs_data[filter_choice]<p_thresh & abs(labs_data$log2FoldChange)>l2FC_thresh,]
  }
  if (label_top_n>0){
    top_res <- DE_res[order(DE_res[filter_choice]),]
    labs_data_top <- top_res[1:label_top_n,]
  }else{
    labs_data_top <- DE_res[FALSE,]
  }

  v_plot <- ggplot(DE_res, aes(x = log2FoldChange, y = -log10(pvalue), color=Significance)) +
    geom_point(size=0.8) +
    geom_hline(yintercept=lowest_sig, linetype="dashed", color = "black")+
    xlab(expression('Log'[2]*'Fold Change')) +
    ylab(expression('-log10(pvalue)')) +
    geom_vline(xintercept = l2FC_thresh, linetype="dashed", color = "black") +
    geom_vline(xintercept = -l2FC_thresh, linetype="dashed", color = "black") +
    ggrepel::geom_label_repel(data = labs_data, mapping = aes(label = gene_id),
                              box.padding = unit(0.35, "lines"),
                              point.padding = unit(0.3, "lines"),
                              force = 1, segment.colour = 'black',show.legend = FALSE,
                              label.size = 1)+

    ggrepel::geom_label_repel(data = labs_data_top, mapping = aes(label = gene_id),
                              box.padding = unit(0.35, "lines"),
                              point.padding = unit(0.3, "lines"),
                              force = 0.5, colour = 'black',show.legend = FALSE) +

    scale_color_manual(breaks = c(paste0("up-reg with ",filter_choice,"<",p_thresh," (",length(sigUp),")"),
                                  paste0("down-reg with ",filter_choice,"<",p_thresh," (",length(sigDown),")"),
                                  paste0("non-significant with ",filter_choice,">",p_thresh," (",length(nonSig),")"),
                                  paste0("low-regulation with ",filter_choice,"<",p_thresh," (",length(lowSig),")")),
                       values=c("#B31B21","#1465AC","darkgray","green")) +
    guides(color = guide_legend(override.aes = list(size = 7)))+
    ggtitle(plot_title) +
    theme_light()+
    theme(text = element_text(size = 10),
          plot.title = element_text(size = 20),
          legend.text = element_text(size = 10))

  return(v_plot)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


#' @title Create an MA plot
#'
#' @description The function creates a 'four quadrant' MA plot, where the FDR and log2FoldChange
#' thresholds dictate the significant up-regulated category, the significant downregulated
#' category, the significant low regulation category and the non-significant category
#'
#' Genes of interest are labeled in a rectangle for visibility and will have the same
#' color as the category in which they are in.
#'
#' The plot is created using ggplot2, to save a plot the ggsave() function is recommended
#' It is also recommended to use the following parameters to save the plot.
#' dpi=300 width=21 height=19 units='cm'
#'
#' @param DE_res The differential expression results to be plotted
#' @param genes_of_interest A vector containing gene names to be labelled in the plot
#' To not label any genes, leave as default or provide an empty vector.
#' @param filter_choice Either padj or pvalue, the choice will be used to filter for
#' significance
#' @param l2FC_thresh The value used as the log2FoldChange threshold
#' @param p_thresh The value used to establish significance from adjusted pvalue
#' @param plot_title The title to be give to the plot
#'
#' @return None
#'
#' @examples
#' TS_object<-create_example_object_for_R()
#' TS_object <- normalize_timeSeries_with_deseq2(time_object=TS_object)
#' #Perform conditional differential gene expression analysis
#' TS_object<-conditional_DE_wrapper(TS_object)
#' DE_res<-slot(TS_object,'DE_results')$conditional$IgM_vs_LPS_TP_1$DE_raw_data
#' ma_plot<-maplot_alt(DE_res = DE_res,filter_choice = 'padj')
#'
#' @import ggplot2
#' @importFrom  ggrepel geom_label_repel
#'
#' @export
#'
maplot_alt <- function(DE_res,genes_of_interest=c(),filter_choice,l2FC_thresh=1,p_thresh=0.05,plot_title="MA plot"){
  #This ensures that the non-significant appear first on the plot (background) while
  #the significant genes are plotted on the foreground

  #Reverse dataframe as to make it descending (highest significance at top of df)
  DE_res<-rev(DE_res)

  DE_res$baseMean_logged <- log2(DE_res$baseMean + 1)
  DE_res$Significance <- NA
  sigUp <- which(DE_res$log2FoldChange>l2FC_thresh & DE_res[filter_choice]<p_thresh)
  DE_res[ sigUp, "Significance"] <- paste0("up-reg | ",filter_choice,"<",p_thresh," (",length(sigUp),")")

  #subset and fill columns based on significant downregulation
  sigDown <- which(DE_res$log2FoldChange<(-l2FC_thresh) & DE_res[filter_choice]<p_thresh)
  DE_res[ sigDown, "Significance"] <- paste0("down-reg | ",filter_choice,"<",p_thresh," (",length(sigDown),")")

  #subset and fill columns based on non-significant upregulation
  nonSig <- which(DE_res[filter_choice]>=p_thresh)
  DE_res[ nonSig, "Significance"] <- paste0("non-significant | ",filter_choice,">",p_thresh," (",length(nonSig),")")

  lowSig <- which(abs(DE_res$log2FoldChange)<=l2FC_thresh & DE_res[filter_choice]<p_thresh)
  DE_res[ lowSig, "Significance"] <- paste0("low-difference | ",filter_choice,"<",p_thresh," (",length(lowSig),")")


  labs_data <- DE_res[DE_res$gene_id %in% genes_of_interest, ]
  ma_plot<-ggplot(DE_res, aes(x=baseMean_logged, y=log2FoldChange,color=Significance))+
    geom_point(size=0.8)+
    geom_hline(yintercept=0, linetype="dashed", color = "red")+
    geom_hline(yintercept=l2FC_thresh, linetype="dashed", color = "black")+
    geom_hline(yintercept=-l2FC_thresh, linetype="dashed", color = "black")+
    xlab(expression('BaseMean')) +
    ylab(expression('log2FoldChange'))+
    ggrepel::geom_label_repel(data = labs_data, mapping = aes(label = gene_id),
                              box.padding = unit(0.35, "lines"),
                              point.padding = unit(0.3, "lines"),
                              force = 1, segment.colour = 'black',show.legend = FALSE,
                              label.size = 0.5) +
    scale_color_manual(breaks = c(paste0("up-reg | ",filter_choice,"<",p_thresh," (",length(sigUp),")"),
                                  paste0("down-reg | ",filter_choice,"<",p_thresh," (",length(sigDown),")"),
                                  paste0("low-difference | ",filter_choice,"<",p_thresh," (",length(lowSig),")"),
                                  paste0("non-significant | ",filter_choice,">",p_thresh," (",length(nonSig),")")),
                       values=c("#B31B21","#1465AC","green","darkgray")) +
    guides(color = guide_legend(override.aes = list(size = 7)))+

    ggtitle(plot_title) +
    theme_light()+
    theme(text = element_text(size = 10),
          plot.title = element_text(size = 20),
          legend.text = element_text(size = 10))

  return(ma_plot)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


#' @title Plot timeseries PCA
#'
#' @description For a DESeq2 analysis, the function which retrieves the variance stabililizing
#' transformation via the \code{vst()} of deseq2.
#' For a limma analysis, the function will calculate the PCs using plotMDS
#'
#' If no DE_type is 'all', the function will plot the PCA for the normalized
#' data.
#'
#' It then plots a labelled PCA plot in png format
#'
#' @param time_object A timeseries object
#' @param exp_name The name of the experiment for which the PCA is plotted
#' @param DE_type Either conditional or temporal for the type of differential experiment
#' being plotted
#' @param show_names boolean indicating if sample names should be put on the pca or not
#'
#' @return the pca_plot
#'
#' @examples
#' TS_object<-create_example_object_for_R()
#' TS_object <- normalize_timeSeries_with_deseq2(time_object=TS_object)
#' #Perform conditional differential gene expression analysis
#' TS_object<-conditional_DE_wrapper(TS_object)
#' TS_object<-temporal_DE_wrapper(TS_object,do_all_combinations=FALSE)
#' TS_pca<-plot_PCA_TS(TS_object,DE_type='all')
#'
#' @import ggplot2
#' @importFrom stats prcomp
#' @importFrom DESeq2 vst plotPCA
#'
#' @export
#'
plot_PCA_TS<-function(time_object,exp_name=NULL,DE_type=NULL,show_names=TRUE){
  samp_dta_full<-exp_sample_data(time_object)

  if(DE_type=='all'){
    DE_res<-NULL
  }else{
    DE_res<-slot(time_object,'DE_results')[[DE_type]][[exp_name]]
  }
  DE_meth<-slot(time_object,'DE_method')
  if(DE_meth=='DESeq2'){
    if (is.null(DE_res)==TRUE){
      dds<-slot(time_object,'DESeq2_obj')
    }else{
      dds<-DE_res[['sub_dds']]
    }
    vsd <- vst(dds, blind=FALSE)
    pca_data <- plotPCA(vsd, intgroup=c('condition'), returnData = TRUE)
    percentVar <- round(100 * attr(pca_data, "percentVar"))

    samp_dta<-samp_dta_full[samp_dta_full$sample %in% pca_data$name,]
    timepoints_used<-samp_dta$timepoint[order(match(samp_dta$sample,pca_data$name))]
    pca_data$timepoint<-factor(timepoints_used,levels=unique(timepoints_used))

  }else if(DE_meth=='limma'){
    if (is.null(DE_res)==TRUE){
      samples_interest<-colnames(exp_matrix(time_object,'norm'))
    }else{
      samples_interest<-colnames(DE_res[['DE_raw_data']])
      samples_interest<-samples_interest[samples_interest %in% samp_dta_full$sample]
    }
    sample_data_used<-samp_dta_full[samp_dta_full$sample %in% samples_interest,]
    sample_data_used<-sample_data_used[order(match(sample_data_used$sample,samples_interest)),]


    matrix_to_use<-exp_matrix(time_object,'norm')[,samples_interest]
    pca_res <- prcomp(t(matrix_to_use))
    axis_labels <- sprintf('%.1f', 1:length(pca_res$sdev), (pca_res$sdev^2 / sum(pca_res$sdev^2))*100)

    percentVar<-c(axis_labels[1],axis_labels[2])

    pca_data<-data.frame(PC1=pca_res$x[,'PC1'],PC2=pca_res$x[,'PC2'],
                         group=sample_data_used$group,
                         timepoint=factor(sample_data_used$timepoint,levels=unique(sample_data_used$timepoint)),
                         name=sample_data_used$sample)

  }
  if(DE_type=='temporal'){
    groups_for_lvl<-strsplit(exp_name,'_vs_')[[1]]
    if(DE_meth=='limma'){
      pca_data$group<-c(paste0('TP_',pca_data$timepoint))
    }
    pca_data$group<-factor(pca_data$group,levels=c(groups_for_lvl[1],groups_for_lvl[2]))
  }else{
    pca_data$group<-factor(pca_data$group,levels=slot(time_object,'group_names'))
  }
  num_shapes<-length(unique(samp_dta_full$timepoint))
  if (is.null(DE_res)==TRUE){
    if(show_names==TRUE){
      pca_plot <- ggplot(pca_data, aes(PC1, PC2, color=group, shape=timepoint,label = name)) +
        geom_label_repel(aes(PC1, PC2, label = name), fontface = 'bold',
                         box.padding = unit(0.35, "lines"),
                         point.padding = unit(0.5, "lines"),
                         segment.color = 'grey50')
    }else{
      pca_plot <- ggplot(pca_data, aes(PC1, PC2, color=group, shape=timepoint))
    }
      pca_plot<-pca_plot+
      scale_shape_manual(values=1:num_shapes) +
      geom_point(size=2, stroke = 2) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance"))
  }else{
    if(show_names==TRUE){
      pca_plot <- ggplot(pca_data, aes(PC1, PC2, color=group, shape=group,label = name)) +
        geom_label_repel(aes(PC1, PC2, label = name), fontface = 'bold',
                         box.padding = unit(0.35, "lines"),
                         point.padding = unit(0.5, "lines"),
                         segment.color = 'grey50')
    }else{
      pca_plot <- ggplot(pca_data, aes(PC1, PC2, color=group, shape=group))
    }
      pca_plot<-pca_plot+
      scale_shape_manual(values=1:num_shapes) +
      geom_point(size=2, stroke = 2) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance"))
  }


  return(pca_plot)
}


# Custom heatmap functions -------------------

#' @title Create heatmap for differential gene expression results
#'
#' @description A wrapper function which calls the necessary functions to create a heatmap
#' illustrating the significant DEGs of a given differential expression type (DE_type)
#'
#' @param time_object A time object
#' @param DE_type The type of heatmap being created (conditional or temporal)
#' @param log_transform Boolean indicating if the results should be log transformed
#' on the heatmap
#' @param plot_file_name The file name to be given to the heatmap once it is saved
#' @param adjust_missing_temp_samples Boolean indicating if missing temporal samples should be
#' compensated for (NA introduced). Currently, there is no reason to have this set to
#' False.
#' @param do_SVGs Boolean indicating if SVG files should be saved for the heatmaps
#' @return none or heatmap plot
#'
#' @examples
#' TS_object<-create_example_object_for_R()
#' TS_object <- normalize_timeSeries_with_deseq2(time_object=TS_object)
#' #Perform conditional differential gene expression analysis
#' TS_object<-conditional_DE_wrapper(TS_object)
#' heat_plot<-custom_heatmap_wrapper(TS_object,DE_type='conditional',log_transform=TRUE,
#'                                   plot_file_name = NULL,
#'                                   adjust_missing_temp_samples=TRUE,do_SVGs=FALSE)
#'
#' @export
#'
custom_heatmap_wrapper<-function(time_object,DE_type,log_transform=TRUE,plot_file_name='custom_DEG_heatmap',adjust_missing_temp_samples=TRUE,do_SVGs=TRUE){

  #Check if sig data exists
  for(name in names(time_object@DE_results[[DE_type]])){
    if(nrow(time_object@DE_results[[DE_type]][[name]][['DE_sig_data']])==0){
      message('Not all comparisons have significant genes, heatmap cannot be produced')
      return(NULL)
    }
  }

  #Create the matrix
  if (DE_type=='conditional'){
    my_heat_list<-create_conditional_heatmap_matrix(time_object)
  }else if (DE_type=='temporal'){
    my_heat_list<-create_temporal_heatmap_matrix(time_object,adjust_for_missing_samples = adjust_missing_temp_samples)
  }

  #Format the data for the plotting function
  temp_list<-prepare_heat_data(my_heat_list,log_transform = log_transform)

  #Associate the results from prepare_heat_data to variables
  my_heat_mat<-temp_list[['heat_matrix']]
  my_region_split<-temp_list[['region_split']]
  my_group_split<-temp_list[['group_split']]
  my_l2fc_vect<-temp_list[['l2fc_vector']]

  #If needed, log transform the l2fc results
  if (log_transform==TRUE){
    my_l2fc_vect<-log_transform_l2fc_vect(my_l2fc_vect)
  }

  DE_meth<-slot(time_object,'DE_method')
  #Call the function to create the heatmap
  if(DE_meth=='limma'){
    legend_val<-'intensity value'
  }else{
    legend_val<-'counts'
  }
  if(is.null(plot_file_name)==TRUE){
    heat_plot<- plot_custom_DE_heatmap(my_heat_mat,my_region_split,my_group_split,
                                       my_l2fc_vect,log_transform = log_transform,
                                       legend_value=legend_val, plot_file_name = plot_file_name,
                                       do_SVG=do_SVGs)
    return(heat_plot)
  }else{
    plot_custom_DE_heatmap(my_heat_mat,my_region_split,my_group_split,my_l2fc_vect,log_transform = log_transform,
                           legend_value=legend_val, plot_file_name = plot_file_name,do_SVG=do_SVGs)
  }

}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @title Create data for conditional heatmap
#'
#' @description This function is intended to be used within the \code{custom_heatmap_wrapper}
#'
#' The function creates a list containing the main data matrix along with three vectors
#' l2fc_vector: A named vector containing the L2FC value for each gene of the matrix
#' group_vector: A vector containing the new names of the groups to represent the group
#' they are in
#' gene_vector: A named vector containing the gene names along with the group that
#' they are in.
#'
#' The vectors are created to account for the possibility of duplicates. For example,
#' a gene may be significant in two timepoint comparisons, and therefore have two
#' different values. The vectors give the genes unique names based on their grouping
#' in order to enable them to have different values if necessary.
#'
#'
#' @param time_object A time series object
#'
#' @return A list containing the main data matrix (count values) and the
#' three vectors described above.
#'
#' @examples
#' TS_object<-create_example_object_for_R()
#' TS_object <- normalize_timeSeries_with_deseq2(time_object=TS_object)
#' #Perform conditional differential gene expression analysis
#' TS_object<-conditional_DE_wrapper(TS_object)
#' heat_dta<-create_conditional_heatmap_matrix(TS_object)
#'
#' @export
#'
create_conditional_heatmap_matrix<-function(time_object){

  sample_dta_full<-exp_sample_data(time_object)
  DE_list<-slot(time_object,'DE_results')$conditional

  DEG_list<-list()
  all_replicates<-c()
  for (i in names(DE_list)){
    if(nrow(DE_list[[i]][['DE_sig_data']])==0){
      next
    }
    DEG_list[[i]]<-DE_list[[i]][['DE_sig_data']]
    replicates_in_exp<-sample_dta_full$replicate[sample_dta_full$sample %in% colnames(DEG_list[[i]])]
    all_replicates<-c(all_replicates,replicates_in_exp)
  }
  rep_template<-unique(all_replicates)

  #Create empty list
  heat_list<-list()

  gene_vector<-c()
  #Create empty df and vector
  heat_df<-as.data.frame(NULL)
  l2fc_vect<-c()

  #Iterate over experiments/timepoints
  for (tp_group in names(DEG_list)){
    # tp_group<-paste0(my_groups,'_TP_',tp)#Create name to look for
    temp_df<-as.data.frame(DEG_list[[tp_group]])#Extract df associated to timepoint

    genes_mod<-temp_df$gene_id#Extract gene names
    genes_mod<-paste0(genes_mod,'_TP',tp_group)#Modify gene names by adding time point

    #Extract and name log2FoldChange
    temp_l2fc<-temp_df$log2FoldChange
    names(temp_l2fc)<-genes_mod

    samp_dta<-sample_dta_full[sample_dta_full$sample %in% colnames(temp_df),]
    if (all(rep_template %in% samp_dta$replicate)==FALSE){
      missing_replicates<-rep_template[!rep_template %in% samp_dta$replicate]
      missing_rep_list<-list()
      for (rep in missing_replicates){
        associated_sample<-sample_dta_full$sample[sample_dta_full$replicate==rep][1]
        missing_rep_list[[associated_sample]]<-rep(NA,nrow(temp_df))
      }
      missing_rep_cols<-as.data.frame(missing_rep_list)
      colnames(missing_rep_cols)=names(missing_rep_list)

      temp_df<-cbind(temp_df,missing_rep_cols)

      replicate_order<-sample_dta_full[sample_dta_full$sample %in% colnames(temp_df),]
      replicate_order<-replicate_order[order(match(replicate_order$sample,colnames(temp_df))),]

      temp_df<-temp_df[,replicate_order$sample]
      colnames(temp_df)=replicate_order$replicate

    }else{
      temp_df<-temp_df[,samp_dta$sample]#Select for experiment related data
      colnames(temp_df)=samp_dta$replicate

    }
    temp_df<-temp_df[,rep_template]#Ensure same order
    row.names(temp_df)=genes_mod
    temp_df<-t(temp_df)#transpose


    #Fill the main dataframe
    if (nrow(heat_df)==0){
      heat_df<-temp_df
    }else{
      heat_df<-cbind(heat_df,temp_df)
    }
    #Concatenate the log2FoldChanges
    l2fc_vect<-c(l2fc_vect,temp_l2fc)
    temp_gene_vector<-rep(paste0(tp_group,' (',length(genes_mod),')'),length(genes_mod))
    gene_vector<-c(gene_vector,temp_gene_vector)
  }

  my_group_df<-sample_dta_full[,c('group','replicate')]

  colnames(my_group_df)=c('labels','replicate')
  my_group_df<-unique(my_group_df)
  row.names(my_group_df)=my_group_df$replicate

  group_vector<-my_group_df[row.names(heat_df),]
  group_vector<-group_vector$labels
  #Fill and return result list
  heat_list[['main_matrix']]<-heat_df
  heat_list[['l2fc_vector']]<-l2fc_vect
  heat_list[['group_vector']]<-group_vector
  heat_list[['gene_vector']]<-gene_vector
  return(heat_list)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @title Create data for temporal heatmap
#'
#' @description This function is intended to be used within the \code{custom_heatmap_wrapper}
#'
#' The function creates a list containing the main data matrix along with three vectors
#' l2fc_vector: A named vector containing the L2FC value for each gene of the matrix
#' group_vector: A vector containing the new names of the groups to represent the group
#' they are in
#' gene_vector: A named vector containing the gene names along with the group that
#' they are in.
#'
#' The vectors are created to account for the possibility of duplicates. For example,
#' a gene may be significant in two timepoint comparisons, and therefore have two
#' different values. The vectors give the genes unique names based on their grouping
#' in order to enable them to have different values if necessary.
#'
#'
#' @param time_object A time series object
#' @param adjust_for_missing_samples Boolean indicating if missing samples should be
#' compensated for (NA introduced). Currently, there is no reason to have this set to
#' False.
#'
#' @return A list containing the main data matrix (count values) and the
#' three vectors described above.
#'
#' @examples
#' TS_object<-create_example_object_for_R()
#' TS_object <- normalize_timeSeries_with_deseq2(time_object=TS_object)
#' #Perform conditional differential gene expression analysis
#' TS_object<-temporal_DE_wrapper(TS_object)
#' heat_dta<-create_temporal_heatmap_matrix(TS_object)
#'
#'
#' @export
#'
create_temporal_heatmap_matrix<-function(time_object,adjust_for_missing_samples=TRUE){
  samp_dta_full<-exp_sample_data(time_object)
  DE_list<-slot(time_object,'DE_results')$temporal
  DEG_list<-list()
  all_replicates<-c()
  for (i in names(DE_list)){
    DEG_list[[i]]<-DE_list[[i]][['DE_sig_data']]
    replicates_in_exp<-samp_dta_full$replicate[samp_dta_full$sample %in% colnames(DEG_list[[i]])]
    all_replicates<-c(all_replicates,replicates_in_exp)
  }
  rep_template<-unique(all_replicates)
  rep_template_ctrl<-paste0(rep_template,'_C')
  rep_template_exp<-paste0(rep_template,'_E')
  rep_template_temporal<-c(rep_template_ctrl,rep_template_exp)

  heat_list<-list()
  gene_vector<-c()
  #Create empty df and vector
  heat_df<-as.data.frame(NULL)
  l2fc_vect<-c()
  for (exp in names(DE_list)){
    ctrl_timepoint<-strsplit(exp,'_')[[1]][5]
    exp_timepoint<-strsplit(exp,'_')[[1]][2]

    temp_df<-DEG_list[[exp]]
    genes_mod<-temp_df$gene_id#Extract gene names
    genes_mod<-paste0(genes_mod,'_',exp)#Modify gene names by adding time point

    #Extract and name log2FoldChange
    temp_l2fc<-temp_df$log2FoldChange
    names(temp_l2fc)<-genes_mod

    sample_data_used<-samp_dta_full[samp_dta_full$sample %in% colnames(temp_df),]

    sample_data_ctrl<-sample_data_used[sample_data_used$timepoint==ctrl_timepoint,]
    sample_data_ctrl$replicate<-paste0(sample_data_ctrl$replicate,'_C')
    sample_data_exp<-sample_data_used[sample_data_used$timepoint==exp_timepoint,]
    sample_data_exp$replicate<-paste0(sample_data_exp$replicate,'_E')

    sample_data_used<-rbind(sample_data_ctrl,sample_data_exp)
    temp_df<-temp_df[,sample_data_used$sample]
    colnames(temp_df)=sample_data_used$replicate

    #Works
    if(adjust_for_missing_samples==TRUE){
      if(all(rep_template_temporal %in% sample_data_used$replicate)==FALSE){
        missing_replicates<-rep_template_temporal[!rep_template_temporal %in% sample_data_used$replicate]
        missing_rep_list<-list()
        for (replicate in missing_replicates){
          missing_rep_list[[replicate]]<-rep(NA,nrow(temp_df))
        }
        missing_rep_cols<-as.data.frame(missing_rep_list)
        colnames(missing_rep_cols)=names(missing_rep_list)

        temp_df<-cbind(temp_df,missing_rep_cols)
      }
      temp_df<-temp_df[,rep_template_temporal]#Ensure format agreement
    }

    row.names(temp_df)=genes_mod
    temp_df<-t(temp_df)#transpose

    #Fill the main dataframe
    if (nrow(heat_df)==0){
      heat_df<-temp_df
    }else{
      heat_df<-cbind(heat_df,temp_df)
    }
    #Concatenate the log2FoldChanges
    l2fc_vect<-c(l2fc_vect,temp_l2fc)
    temp_gene_vector<-rep(paste0(exp,' (',length(genes_mod),')'),length(genes_mod))
    gene_vector<-c(gene_vector,temp_gene_vector)

  }


  group_vector<-c()
  for (group in row.names(heat_df)){
    if (endsWith(group,'C')==TRUE){
      group_vector<-c(group_vector,'Control')
    }else{
      group_vector<-c(group_vector,'Experiment')
    }
  }

  #Fill and return result list
  heat_list[['main_matrix']]<-heat_df
  heat_list[['l2fc_vector']]<-l2fc_vect
  heat_list[['group_vector']]<-group_vector
  heat_list[['gene_vector']]<-gene_vector

  return(heat_list)
}
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @title Process data for DE heatmap
#'
#' @description Function which adjusts and renames certain elements of the heat_list
#' produced by \code{create_heatmap_matrix}
#'
#' It performs a log_transformation if necessary and converts the vectors into
#' factors to be used for annotation purposes within the heatmap
#'
#' This function is intended to be used within the \code{cusom_heatmap_wrapper}
#'
#' @param matrix_list The heat_list produced by \code{create_heatmap_matrix}
#' @param log_transform If a log transformation of the data is necessary
#'
#' @return An updated version of the inputted matrix_list
#'
#' @examples
#' TS_object<-create_example_object_for_R()
#' TS_object <- normalize_timeSeries_with_deseq2(time_object=TS_object)
#' #Perform conditional differential gene expression analysis
#' TS_object<-conditional_DE_wrapper(TS_object)
#' heat_dta<-create_conditional_heatmap_matrix(TS_object)
#' heat_dta<-prepare_heat_data(heat_dta,log_transform=TRUE)
#'
#' @export
#'
prepare_heat_data <- function(matrix_list,log_transform){
  heat_matrix<-matrix_list[['main_matrix']]

  if (log_transform==TRUE){
    heat_matrix<-log10(as.matrix(heat_matrix+1))
  }

  group_split<-matrix_list[['group_vector']]

  numbered_regions<-matrix_list[['gene_vector']]

  #Sets up the order of appearance
  region_split <- factor(numbered_regions, levels=unique(numbered_regions))
  group_split <- factor(group_split, levels=unique(group_split))

  l2fc_vector<-matrix_list[['l2fc_vector']]
  names(l2fc_vector)=region_split
  return_list<-list('heat_matrix'=heat_matrix,'region_split'=region_split,
                    'group_split'=group_split,'l2fc_vector'=l2fc_vector)
  return(return_list)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @title Log transform L2FC/FC values
#'
#' @description The function first identifies the location of negative values. It then performs
#' the log transformation on the absolute values of the vector.
#' Lastly it returns originally negative values to negative values
#' This is done as log transformation does not work on negative values, therefore a
#' work around had to be done.
#'
#' This function is intended to be used within \code{custom_heatmap_wrapper}
#'
#' @param l2fc_vector The log2FoldChange vector to be log transformed
#'
#' @return The log transformed vector
#'
#' @examples
#' TS_object<-create_example_object_for_R()
#' TS_object <- normalize_timeSeries_with_deseq2(time_object=TS_object)
#' #Perform conditional differential gene expression analysis
#' TS_object<-conditional_DE_wrapper(TS_object)
#' heat_dta<-create_conditional_heatmap_matrix(TS_object)
#' heat_dta<-prepare_heat_data(heat_dta,log_transform=TRUE)
#'
#' log_l2fc<-log_transform_l2fc_vect(heat_dta[['l2fc_vector']])
#'
#' @export
#'
log_transform_l2fc_vect <-function(l2fc_vector){
  neg_idx<-c()
  for (idx in 1:length(l2fc_vector)){
    if (l2fc_vector[idx]<0){
      neg_idx<-c(neg_idx,idx)
    }
  }


  new_l2fc<-log10(abs(l2fc_vector))
  for (idx in 1:length(new_l2fc)){
    if (idx %in% neg_idx){
      new_l2fc[idx]<- -new_l2fc[idx]
    }
  }

  return(new_l2fc)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @title Plots DE heatmaps
#'
#' @description The function is intended to be called within the \code{custom_heatmap_wrapper}
#'
#' The function creates a segmented heatmap where rows are samples/patients and columns
#' are genes. The heatmap is column segmented based on the number of differential expression
#' analyses contained within it. A log2foldChange histogram is attached to the bottom
#' of the heatmap to indicate significance of the genes being visualized.
#'
#' The column/segment legend is given at the bottom of the heatmap with the number
#' of genes in each group being also provided.
#'
#' The gradient legend of the values within the heatmap is given on the right hand side.
#'
#' The heatmap is saved in both png and svg format.
#'
#'
#' @param heat_mat The main count matrix of the heatmap
#' @param col_split The named vector used to split the columns into it's segments
#' @param row_splits The named vector used to split the rows into the two groups
#' @param l2fc_col The log2foldchange vector used to create the histogram
#' @param log_transform boolean indicating if the results were log transformed or not
#' @param legend_value string indicating the type of value used
#' (counts for RNAseq, intensity value for microarray)
#' @param plot_file_name The name given to the saved heatmap
#' @param custom_width The width of the heatmap
#' @param custom_height The height of the heatmap
#' @param do_SVG Boolean to check if the SVG file should be created or not
#'
#' @return none or heatmap plot
#'
#' @examples
#' TS_object<-create_example_object_for_R()
#' TS_object <- normalize_timeSeries_with_deseq2(time_object=TS_object)
#' #Perform conditional differential gene expression analysis
#' TS_object<-conditional_DE_wrapper(TS_object)
#' heat_dta<-create_conditional_heatmap_matrix(TS_object)
#' heat_dta<-prepare_heat_data(heat_dta,log_transform=TRUE)
#'
#' log_l2fc<-log_transform_l2fc_vect(heat_dta[['l2fc_vector']])
#'
#' heat_plot<- plot_custom_DE_heatmap(heat_dta[['heat_matrix']],heat_dta[['region_split']],
#'                                    heat_dta[['group_split']],log_l2fc,log_transform = TRUE,
#'                                    legend_value='counts', plot_file_name = NULL,
#'                                    do_SVG=FALSE)
#'
#' @importFrom ComplexHeatmap HeatmapAnnotation anno_barplot rowAnnotation anno_block Heatmap draw Legend
#' @importFrom grid gpar
#' @export
#'
plot_custom_DE_heatmap <-function(heat_mat,col_split,row_splits,l2fc_col, log_transform,
                                  legend_value='counts', plot_file_name='custom_heatmap',
                                  custom_width=15,custom_height=5,do_SVG=TRUE){

  if(log_transform==TRUE){
    if(legend_value=='intensity value'){
      bottom_histo = HeatmapAnnotation('log10(FC)' = anno_barplot(l2fc_col))
    }else{
      bottom_histo = HeatmapAnnotation('log10(l2FC)' = anno_barplot(l2fc_col))
    }
    count_legend<-paste0('log10(',legend_value,')')
  }else{
    if(legend_value=='intensity value'){
      bottom_histo = HeatmapAnnotation('FC' = anno_barplot(l2fc_col))
    }else{
      bottom_histo = HeatmapAnnotation('l2FC' = anno_barplot(l2fc_col))
    }
    count_legend<-legend_value
  }

  fill_set_regions=gpar(fill=2:(length(unique(col_split))+1))
  top_annot = HeatmapAnnotation(foo = anno_block(gp = fill_set_regions, labels = rep(NULL,length(unique(col_split)))))

  fill_set_groups=gpar(fill = c("#998ec3", "#f1a340"))
  left_annot = rowAnnotation(foo = anno_block(gp = fill_set_groups, labels = unique(row_splits)))


  legend_label<-unique(col_split)
  #Calculate ideal distance between legend points in mm, the max distance will be used
  gap_vect<-nchar(levels(legend_label))*1.9

  #Create legend
  lgd = Legend(labels = legend_label, title = "Experiments", legend_gp = fill_set_regions,nrow = 1,
               title_position = "leftcenter",gap = unit(max(gap_vect), "mm"))

  #Check if clustering of rows and cols is possible
  problematic_combs<-identify_problematic_combs(heat_mat)
  if (length(problematic_combs$row)>0){
    clust_rows<-FALSE
  }else{
    clust_rows<-TRUE
  }

  if (length(problematic_combs$column)>0){
    clust_cols<-FALSE
  }else{
    clust_cols<-TRUE
  }

  if(is.null(plot_file_name)==TRUE){
    pdf(NULL)
    heat_plot<-draw(Heatmap(heat_mat,name=count_legend,cluster_rows = clust_rows,cluster_columns = clust_cols,
                              show_column_names = FALSE,show_row_names = FALSE,row_names_side='left',
                              row_split = row_splits,column_split=col_split,
                              border=TRUE,na_col = 'gray',column_title = NULL,row_title = NULL,
                              top_annotation = top_annot,left_annotation=left_annot,
                              bottom_annotation = bottom_histo,
                              cluster_column_slices=FALSE,cluster_row_slices = FALSE),
                      annotation_legend_list = lgd,
                      annotation_legend_side = 'bottom'
    )
    trash<-capture.output(dev.off())#Capture output to prevent print
    return(heat_plot)
  }

  if(do_SVG==TRUE){
    save_name_svg<-paste0(plot_file_name,'.svg')
    svg(save_name_svg,width=custom_width,height=custom_height)
    draw(Heatmap(heat_mat,name=count_legend,cluster_rows = clust_rows,cluster_columns = clust_cols,
                 show_column_names = FALSE,show_row_names = FALSE,row_names_side='left',
                 row_split = row_splits,column_split=col_split,
                 border=TRUE,na_col = 'gray',column_title = NULL,row_title = NULL,
                 top_annotation = top_annot,left_annotation=left_annot,
                 bottom_annotation = bottom_histo,
                 cluster_column_slices=FALSE,cluster_row_slices = FALSE),
         annotation_legend_list = lgd,
         annotation_legend_side = 'bottom'
    )
    dev.off()
  }

  save_name_png<-paste0(plot_file_name,'.png')

  #width and height *96 to convert inches to pixels
  png(save_name_png,width=custom_width*96,height=custom_height*96)
  draw(Heatmap(heat_mat,name=count_legend,cluster_rows = clust_rows,cluster_columns = clust_cols,
               show_column_names = FALSE,show_row_names = FALSE,row_names_side='left',
               row_split = row_splits,column_split=col_split,
               border=TRUE,na_col = 'gray',column_title = NULL,row_title = NULL,
               top_annotation = top_annot,left_annotation=left_annot,
               bottom_annotation = bottom_histo,
               cluster_column_slices=FALSE,cluster_row_slices = FALSE),
       annotation_legend_list = lgd,
       annotation_legend_side = 'bottom'
  )
  dev.off()
}



# PART heatmap function -------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @title Create top annotation for PART heatmap
#'
#' @description This function is intended to be used within \code{PART_heat_map}
#'
#' The top annotation shows annotation blocks for the groups as well as annotation
#' gradients for the time points.
#' Two sets of top annotations are created, one labelled and one unlabeled.
#' In addition, the column splits are created based on the replicates
#'
#'
#' @param object A time series object
#'
#' @return A list containing the two types of top annotaiton and the named vector used to create the column splits
#'
#' @examples
#' TS_object<-create_example_object_for_R()
#' TS_object <- normalize_timeSeries_with_deseq2(time_object=TS_object)
#' #Perform conditional differential gene expression analysis
#' TS_object<-conditional_DE_wrapper(TS_object)
#'
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
#' top_annot<-prepare_top_annotation_PART_heat(TS_object)
#'
#' @importFrom ComplexHeatmap HeatmapAnnotation
#' @importFrom grid gpar
#' @import RColorBrewer
#' @import grDevices
#'
#' @export
#'
prepare_top_annotation_PART_heat<-function(object){

  main_matrix<-slot(object,'PART_results')$part_matrix
  samp_dta_full<-exp_sample_data(object)
  group_cols<-slot(object,'group_colors')
  samp_data<-samp_dta_full[order(match(samp_dta_full$sample,colnames(main_matrix))),]

  found_timepoints<-samp_data$timepoint
  found_replicates<-samp_data$replicate
  group_order<-unique(samp_data$group)

  col_split<-factor(unname(found_replicates),levels=unique(unname(found_replicates)))
  num_cols<-table(sapply(strsplit(levels(col_split),"_[0-9]"), `[`, 1))
  #Needs to be adatapted using the object graphical parameters
  fill_set_groups=gpar(fill=c(rep(group_cols[[group_order[1]]],unname(num_cols[1])),
                              rep(group_cols[[group_order[2]]],unname(num_cols[2]))
  ))

  my_cols<-c('#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026')
  # If possible, use custom colors, otherwise use colorbrewer
  if (length(unique(found_timepoints))>length(my_cols)){
    num_tp <- length(unique(found_timepoints))
    cols <- brewer.pal(8, name = "Set1")[seq_len(min(8, num_tp))]
    time_cols <- colorRampPalette(colors = cols)(num_tp)
    names(time_cols) <- unique(found_timepoints)
  }else{
    time_cols<-my_cols[1:length(unique(found_timepoints))]
    names(time_cols)<-unique(found_timepoints)
  }


  top_annot_labels = HeatmapAnnotation(foo = anno_block(gp = fill_set_groups, labels = unique(col_split)),
                                       timepoints=found_timepoints,
                                       col = list(timepoints = time_cols),
                                       show_annotation_name=FALSE,
                                       annotation_legend_param = list(title = "timepoints", at = unique(found_timepoints),
                                                                      labels = unique(found_timepoints)))

  top_annot_no_labels = HeatmapAnnotation(foo = anno_block(gp = fill_set_groups, labels = rep(NULL,length(unique(col_split)))),
                                          timepoints=found_timepoints,
                                          col = list(timepoints = time_cols),
                                          show_annotation_name=FALSE,
                                          annotation_legend_param = list(title = "timepoints", at = unique(found_timepoints),
                                                                         labels = unique(found_timepoints)))

  return(list(top_annot_labels,top_annot_no_labels,col_split))
}

#' @title Create and plot PART heatmap
#'
#' @description The heatmap shows selected genes as rows and replicates of the time series experiment
#' as columns.
#' Gene clusters are identified by colored row annotation. Replicates are split
#' by their grouping (conditional) and ordered (within each other) via time points
#' The legend for groupings, clusters, timepoints, and heatmap values is given on
#' the right hand side of the heatmap.
#'
#' the heatmap is saved twice, once
#' in png format, and the other in svg format
#'
#'
#' @param object A time series object
#' @param heat_name The file name given to the saved heatmap
#'
#' @importFrom ComplexHeatmap rowAnnotation Heatmap draw
#' @importFrom grid gpar
#'
#' @return none or PART heatmap if heat_name is set to NULL
#'
#' @examples
#' TS_object<-create_example_object_for_R()
#' TS_object <- normalize_timeSeries_with_deseq2(time_object=TS_object)
#' #Perform conditional differential gene expression analysis
#' TS_object<-conditional_DE_wrapper(TS_object)
#'
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
#' #Heatmap will be saved to main directory
#' PART_heat<-PART_heat_map(TS_object,NULL) #Create a summary heatmap
#'
#' @export
#'
PART_heat_map<-function(object, heat_name='custom_heat_map'){

  #Cluster illustration stored as rowAnnotation
  PART_res<-slot(object,'PART_results')
  row_annot <- rowAnnotation(gene_cluster = PART_res$part_data$gene_cluster,
                             col = list(gene_cluster=PART_res$cluster_info[['colored_clust_rows']]),
                             show_annotation_name=FALSE,
                             annotation_legend_param = list(title = "clusters", at = unique(PART_res$part_data$gene_cluster),
                                                            labels = unique(PART_res$cluster_map$cluster)))

  #Create top annotations
  top_annot_results<-prepare_top_annotation_PART_heat(object)
  top_annot_labels<-top_annot_results[[1]]
  top_annot_no_labels<-top_annot_results[[2]]
  col_split<-top_annot_results[[3]]

  #Prepare a 'gap_vect' to split the different elements of the heatmap
  col_split_vect<-levels(col_split)
  gap_vect<-c()
  for (idx in 1:length(col_split_vect)){
    if (idx!=length(col_split_vect)){
      val1<-unlist(strsplit(as.character(col_split_vect[idx]),'_'))[1]
      val2<-unlist(strsplit(as.character(col_split_vect[idx+1]),'_'))[1]
      if (val1!=val2){
        gap_vect<-c(gap_vect,4)
      }else{
        gap_vect<-c(gap_vect,0.5)
      }
    }
  }

  group_cols<-slot(object,'group_colors')
  #Create the legend for the groups
  lgd = Legend(labels = names(group_cols), title = "groups",
               legend_gp = gpar(fill = unname(group_cols)))


  #Extract matrix for plotting
  sorted_matrix<-as.matrix(PART_res$part_data[,3:ncol(PART_res$part_data)])

  #If save location is NULL, return the plot instead of saving
  if(is.null(heat_name)==TRUE){
    pdf(NULL)
    PART_plot<- draw(
      Heatmap(
        sorted_matrix, name = "Z-score", cluster_columns = FALSE,
        cluster_rows=FALSE,#PART_res$cluster_info[['clustered_rows']],
        show_column_dend = TRUE,show_row_dend = FALSE,
        row_names_gp = gpar(fontsize = 8), left_annotation = row_annot,
        row_order = row.names(PART_res$part_data),
        show_row_names = FALSE,top_annotation = top_annot_no_labels,column_split = col_split,cluster_column_slices = TRUE,
        column_gap = unit(gap_vect, "mm"),show_column_names = FALSE,border=FALSE,column_title = NULL),
      annotation_legend_list = lgd
    )
    dev.off()
    return(PART_plot)
  }else{
    PART_save_data<-paste0(heat_name,'_data.csv')
    PART_save_cmap<-paste0(heat_name,'_cmap.csv')
    write.csv(PART_res$part_matrix,PART_save_data)
    write.csv(PART_res$cluster_map,PART_save_cmap)
  }

  #Plot the heatmap
  svg(paste0(heat_name,"_with_names.svg"),height=30,width=20)
  draw(
    Heatmap(
      sorted_matrix, name = "Z-score", cluster_columns = FALSE,
      cluster_rows=FALSE,#PART_res$cluster_info[['clustered_rows']],
      show_column_dend = TRUE,show_row_dend = FALSE,
      row_names_gp = gpar(fontsize = 8), left_annotation = row_annot,
      row_order = row.names(PART_res$part_data),
      show_row_names = TRUE,top_annotation = top_annot_labels,column_split = col_split,cluster_column_slices = TRUE,
      column_gap = unit(gap_vect, "mm"),show_column_names = FALSE,border=FALSE,column_title = NULL),
    annotation_legend_list = lgd
  )
  trash<-capture.output(dev.off())#Capture output to prevent print


  svg(paste0(heat_name,".svg"))
  draw(
    Heatmap(
      sorted_matrix, name = "Z-score", cluster_columns = FALSE,
      cluster_rows=FALSE,#PART_res$cluster_info[['clustered_rows']],
      show_column_dend = TRUE,show_row_dend = FALSE,
      row_names_gp = gpar(fontsize = 8), left_annotation = row_annot,
      row_order = row.names(PART_res$part_data),
      show_row_names = FALSE,top_annotation = top_annot_no_labels,column_split = col_split,cluster_column_slices = TRUE,
      column_gap = unit(gap_vect, "mm"),show_column_names = FALSE,border=FALSE,column_title = NULL),
    annotation_legend_list = lgd
  )
  trash<-capture.output(dev.off())#Capture output to prevent print
}


# cluster trajectory functions -------------------

#' @title Calcualte trajectory of clusters
#'
#' @description Function which calculates the trajectory data for each cluster
#' Trajectory data is defined as the impact a gene has on time
#' This is obtained by performing a scale feature sum, where each
#' gene is 'equalized' by dividing the value of the gene (at each sample)
#' by the rowSum of the gene (the addition of the gene's values across all samples)
#'
#' @param object A timeseries object
#' @param custom_cmap A custom cluster map to be used instead of the cmap contained within the object
#' @param scale_feat If the genes should be scaled/transformed using scale feature sum
#'
#' @return A dataframe containing the transformed or non-transformed gene values for
#' each cluster
#'
#' @examples
#' TS_object<-create_example_object_for_R()
#' TS_object <- normalize_timeSeries_with_deseq2(time_object=TS_object)
#' #Perform conditional differential gene expression analysis
#' TS_object<-conditional_DE_wrapper(TS_object)
#'
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
#' ts_data<-calculate_cluster_traj_data(TS_object,scale_feat=TRUE) #Calculate scaled gene values for genes of clusters
#'
#' @importFrom reshape2 melt
#' @importFrom stringi stri_reverse
#' @importFrom stringr str_split_fixed
#'
#' @export
#'
calculate_cluster_traj_data<-function(object,custom_cmap=NULL,scale_feat=TRUE){
  PART_res<-slot(object,'PART_results')
  samp_dta_full<-exp_sample_data(object)
  if (is.null(custom_cmap)==TRUE){
    my_cmap<-PART_res$cluster_map
  }else{
    my_cmap<-custom_cmap
  }

  norm_mat<-exp_matrix(object,'norm')[row.names(my_cmap),]

  ts_df<-data.frame(gene_id=NULL,group=NULL,timepoint=NULL,mean_reads=NULL)
  df_list<-list()
  for (group in unique(samp_dta_full$group)){
    for (tp in unique(samp_dta_full$timepoint)){
      samples_to_merge<-samp_dta_full$sample[samp_dta_full$group==group & samp_dta_full$timepoint==tp]
      if(length(samples_to_merge)>1){
        val_vect<-rowMeans(norm_mat[,samples_to_merge])
      }else{
        val_vect<-norm_mat[,samples_to_merge,drop=FALSE]
      }
      new_col_name<-paste0(group,'_',tp)
      df_list[[new_col_name]]<-val_vect
    }
  }
  if(scale_feat==TRUE){
    ts_df<-sweep(data.frame(df_list),1,rowSums(data.frame(df_list)),'/')
  }else{
    ts_df<-(data.frame(df_list))
  }
  ts_df$gene_id<-row.names(ts_df)
  ts_df <- melt(ts_df,id='gene_id')
  for_merge<-samp_dta_full[,c('sample','group','timepoint')]

  # colnames(ts_df)=c('gene_id','sample','trans_mean')
  # print(ts_df)
  # print(for_merge)
  # ts_df<-merge(ts_df,for_merge,by='sample')
  #
  # ts_df<-ts_df[,c('gene_id','group','timepoint','trans_mean')]

  # the column is reversed before splitting to ensure that the '_' being split
  # on is the one before the timepoint
  ts_df<-cbind(ts_df,str_split_fixed(stri_reverse(ts_df$variable), "_", 2))

  colnames(ts_df)=c('gene_id','replicate','trans_mean','timepoint','group')

  #Re-reverse the groups and timepoints to ensure proper order following split
  ts_df$group<-stri_reverse(ts_df$group)
  ts_df$timepoint<-stri_reverse(ts_df$timepoint)

  ts_df<-ts_df[,c('gene_id','group','timepoint','trans_mean')]
  ts_df$timepoint<-as.numeric(ts_df$timepoint)
  #Add cluster association
  my_cmap$gene_id<-row.names(my_cmap)
  my_cmap<-my_cmap[,c('gene_id','cluster')]
  clust_df<-as.data.frame(table(my_cmap$cluster))
  colnames(clust_df)=c('cluster','nGenes')
  clust_df$labels<-paste0(clust_df$cluster,' - ',clust_df$nGenes,' genes |')

  ts_df<-merge(ts_df,my_cmap,by='gene_id')
  ts_df<-merge(ts_df,clust_df,by='cluster')
  # View(ts_df)
  ts_df$labels<-paste0(ts_df$labels,' ',ts_df$group)
  #Order the dataframe
  ts_df<-ts_df[order(ts_df$cluster,ts_df$group,ts_df$gene_id,ts_df$timepoint),]

  return(ts_df)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


#' @title Calculate mean trajectory of clusters
#'
#' @description Function which calculates the mean trajectory of each cluster
#'
#' @param clust_traj_dta gene trajectory data for each gene, obtained by \code{calculate_cluster_traj_data}
#'
#' @return A dataframe containing the mean cluster trajectories
#'
#' @examples
#' TS_object<-create_example_object_for_R()
#' TS_object <- normalize_timeSeries_with_deseq2(time_object=TS_object)
#' #Perform conditional differential gene expression analysis
#' TS_object<-conditional_DE_wrapper(TS_object)
#'
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
#' ts_data<-calculate_cluster_traj_data(TS_object,scale_feat=TRUE) #Calculate scaled gene values for genes of clusters
#' mean_ts_data<-calculate_mean_cluster_traj(ts_data) #Calculate the mean scaled values for each cluster
#'
#' @export
#'
calculate_mean_cluster_traj<-function(clust_traj_dta){
  mean_clust_df<-data.frame(NULL)
  for (group in unique(clust_traj_dta$group)){
    for (clust in unique(clust_traj_dta$cluster)){
      for (tp in unique(clust_traj_dta$timepoint)){
        sub_df<-clust_traj_dta[clust_traj_dta$group==group &
                                 clust_traj_dta$cluster==clust &
                                 clust_traj_dta$timepoint==tp,]

        mean_val<-mean(sub_df$trans_mean)
        temp_df<-data.frame(group=group,cluster=clust,timepoint=tp,
                            trans_mean=mean_val,labels=unique(sub_df$labels))
        if(nrow(mean_clust_df)==0){
          mean_clust_df<-temp_df
        }else{
          mean_clust_df<-rbind(mean_clust_df,temp_df)
        }
      }
    }
  }
  mean_clust_df<-mean_clust_df[order(mean_clust_df$cluster,mean_clust_df$group,mean_clust_df$timepoint),]
  return(mean_clust_df)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


#' @title Plot cluster trajectories
#'
#' @description Function which creates a grid of cluster trajectories. Each cluster is split
#' into two subplots, one for the control and the other for the experiment.
#'
#' Individual gene trajectories are plotted for each cluster along with a large
#' gray line for the mean cluster trajectories.
#'
#' @param object A timeseries object
#' @param ts_data The trajectory data for all clusters being calculated
#' The data is calculated/obtained from \code{calculate_cluster_traj_data} function
#' @param ts_mean_data The trajectory data for all clusters being calculated
#' The data is calculated/obtained from \code{calculate_cluster_traj_data} function
#' @param num_col Integer stating the number of columns for the plots.
#' @param rem_legend_axis Boolean indicating if the legend and axis titles should be removed
#' @param log_TP Boolean indicating if timepoints should be log transformed
#' @param title_text_size Integer indicating what the font size of the titles should be in the facets
#'
#' @return A ggplot2 object for the cluster trajectory plot performed
#'
#' @examples
#' TS_object<-create_example_object_for_R()
#' TS_object <- normalize_timeSeries_with_deseq2(time_object=TS_object)
#' #Perform conditional differential gene expression analysis
#' TS_object<-conditional_DE_wrapper(TS_object)
#'
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
#' ts_data<-calculate_cluster_traj_data(TS_object,scale_feat=TRUE) #Calculate scaled gene values for genes of clusters
#' mean_ts_data<-calculate_mean_cluster_traj(ts_data) #Calculate the mean scaled values for each cluster
#' clust_traj<-plot_cluster_traj(TS_object,ts_data,mean_ts_data)
#'
#' @import ggplot2
#'
#' @export
#'
plot_cluster_traj<-function(object,ts_data,ts_mean_data,num_col=4,rem_legend_axis=FALSE,log_TP=FALSE,title_text_size=14){

  if(log_TP==TRUE){
    ts_data$log10_timepoint<-log10(ts_data$timepoint)
    ts_data$log10_timepoint[ts_data$log10_timepoint=='-Inf']<-0
    ts_mean_data$log10_timepoint<-log10(ts_mean_data$timepoint)
    ts_mean_data$log10_timepoint[ts_mean_data$log10_timepoint=='-Inf']<-0
  }


  if('log10_timepoint' %in% colnames(ts_data)){
    plt <- ggplot(ts_data, aes(y = trans_mean , x = `log10_timepoint`, color = group))
  }else{
    plt <- ggplot(ts_data, aes(y = trans_mean , x = timepoint, color = group))
  }
  plt <- plt + scale_color_manual(values=slot(object,'group_colors')) +
    geom_line(aes(group = gene_id), alpha = 0.4) +
    geom_point() +
    geom_line(
      data = ts_mean_data, lwd = 1.5, color = "grey50",
      aes(group = group)
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    # scale_y_continuous(expand = c(0, 0))+
    ylab('scaled expression') +
    facet_wrap(~labels, scales = 'free_x', ncol = num_col)

  if(rem_legend_axis==TRUE){
    plt<-plt + theme(legend.position = "none") +
      theme(axis.title.x = element_blank()) +
      theme(axis.title.y = element_blank())
  }
  plt <- plt + theme(strip.text.x = element_text(size = title_text_size))

  return(plt)
}

# gene trajectory functions -------------------

#' @title Calculate trajectory of a single gene
#'
#' @description Function takes in a target gene and will calculate the genes trajectory
#' over the time points within the time object
#'
#' @param time_object A timeseries object
#' @param target_gene String indicating the gene
#' @param log_timepoint Boolean indicating if timepoints should be log transformed
#'
#' @return dataframe with the trajectory data for the requested gene
#'
#' @examples
#' TS_object<-create_example_object_for_R()
#' TS_object <- normalize_timeSeries_with_deseq2(time_object=TS_object)
#' #Perform conditional differential gene expression analysis
#' TS_object<-conditional_DE_wrapper(TS_object)
#' aicda_traj_dta<-calculate_gene_traj_data(TS_object,'AICDA')
#'
#' @importFrom reshape2 melt
#'
#' @export
calculate_gene_traj_data<-function(time_object,target_gene,log_timepoint=FALSE){
  #Create a list of genes of interest
  my_dta<-exp_matrix(time_object,'norm')[target_gene,]

  my_dta<-melt(my_dta,id.vars = NULL)
  my_dta$sample<-row.names(my_dta)
  if(slot(time_object,'DE_method')=='limma'){
    colnames(my_dta)=c('reads','sample')
  }else{
    colnames(my_dta)=c('sample','reads')
  }
  sample_data<-exp_sample_data(time_object)
  my_dta<-merge(my_dta,sample_data,by='sample')
  mean_data_list<-list()
  for (group in unique(my_dta$group)){
    for (tp in unique(my_dta$timepoint)){
      new_val=mean(my_dta$reads[my_dta$group==group & my_dta$timepoint==tp])
      mean_data_list[['reads']]<-c(mean_data_list[['reads']],new_val)
      mean_data_list[['timepoint']]<-c(mean_data_list[['timepoint']],tp)
      mean_data_list[['label']]<-c(mean_data_list[['label']],target_gene)
      mean_data_list[['group']]<-c(mean_data_list[['group']],group)
    }
  }
  mean_data<-as.data.frame(mean_data_list)
  if(log_timepoint==TRUE){
    mean_data$log10_timepoint<-log10(mean_data$timepoint)
    mean_data$log10_timepoint[mean_data$log10_timepoint=='-Inf']<-0
  }
  return(mean_data)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


#' @title Plot gene trajectory
#'
#' @description The function plots the gene's trajectory across all timepoints for all the
#' groups contained in the dataframe. It performs a facet split based on the 'label'
#' column of the dataframe
#'
#' @param mean_data dataframe containing mean reads per timepoint
#' @param color_vector a vector with color IDs where the names match the groups names
#' (by default it is Activated and Control)
#'
#' @return the ggplot2 object for the plot
#'
#' @examples
#' TS_object<-create_example_object_for_R()
#' TS_object <- normalize_timeSeries_with_deseq2(time_object=TS_object)
#' #Perform conditional differential gene expression analysis
#' TS_object<-conditional_DE_wrapper(TS_object)
#' aicda_traj_dta<-calculate_gene_traj_data(TS_object,'AICDA')
#' group_cols<-slot(TS_object,'group_colors')
#' aicda_plot<-plot_single_gene_traj(aicda_traj_dta,group_cols)
#'
#' @import ggplot2
#'
#' @export
#'
plot_single_gene_traj<-function(mean_data,color_vector=NULL){
  if(is.null(color_vector)==TRUE){
    color_vector<-c("#e31a1c","#1f78b4")
    names(color_vector)=c('Activated','Control')
  }
  if('log10_timepoint' %in% colnames(mean_data)){
    plt <- ggplot(mean_data,aes(x = log10_timepoint, y = reads, color = group))+
      geom_smooth(aes(x = log10_timepoint ), lwd = 1.5,se=FALSE,formula = y ~ x, method = "loess")
  }else{
    plt <- ggplot(mean_data,aes(x = timepoint, y = reads, color = group))+
      geom_smooth(aes(x = timepoint ), lwd = 1.5,se=FALSE,formula = y ~ x, method = "loess")
  }
  plt <- plt +
    scale_color_manual(values=color_vector) +
    geom_point(size = 3) +
    # scale_x_continuous(expand = c(0, 0)) +
    facet_wrap(~ label, scales = 'free_x', ncol = 1)

  if(length(unique(mean_data$timepoint))<3){
    plt<-plt+geom_line(aes(group = group))
  }

  return(plt)
}



# data creation functions -------------------

#' @title Create and save DE results
#'
#' @description Function which saves the raw and significant version of the differential
#' gene expression data for a specific differential gene expression experiment
#' within a timeseries object
#'
#' @param object A timeseries object
#' @param DE_type Character for the type of differential gene expression containing
#' the experiment of interest (conditional or temporal)
#' @param exp_name Character stating the name of the differential gene expression
#' of interest
#' @param save_location Character the location where the csv files will be saved.
#'
#' @return None
#'
#' @examples
#' TS_object<-create_example_object_for_R()
#' TS_object <- normalize_timeSeries_with_deseq2(time_object=TS_object)
#' #Perform conditional differential gene expression analysis
#' TS_object<-conditional_DE_wrapper(TS_object)
#' #Below function will save results to main directory
#' my_res<-create_DE_data_results(TS_object,DE_type='conditional',exp_name='IgM_vs_LPS_TP_1',save_location=NULL)
#'
#' @export
#'
create_DE_data_results<-function(object,DE_type,exp_name,save_location){
  DE_res<-slot(object,'DE_results')
  DE_raw<-DE_res[[DE_type]][[exp_name]][['DE_raw_data']]
  DE_sig<-DE_res[[DE_type]][[exp_name]][['DE_sig_data']]

  if(is.null(save_location)==TRUE){
    return_list<-list('DE_raw'=DE_raw,'DE_sig'=DE_sig)
  }else{
    write.csv(DE_raw,paste0(save_location,'DE_raw_data.csv'))
    write.csv(DE_sig,paste0(save_location,'DE_sig_data.csv'))
  }
}


#' @title Creates tables and plots for genes of interest
#'
#' @description Function which retrieves, processes, and saves information regarding
#' genes which have been labelled as interesting by the user. The function will
#' retrieve the genes differential expression results on both the conditional and temporal
#' axes, then save these results in csv format. One csv file per gene of interest.
#'
#' The function will also plot and save (as png) the trajectory of each gene
#'
#' @param object A TimeSEries object
#' @param genes_of_interest A vector of strings indicating the genes of interest
#' @param save_location The directory in which the results will be saved
#' @param DE_type The differential gene expression axes for which the results will be retrieved
#' @param log_tp Boolean indicating if timepoints should be log10 transformed for the plot
#'
#' @return None
#'
#' @import ggplot2
#'
#' @export
#'
create_tables_genes_of_interest_DE<-function(object,genes_of_interest,save_location,DE_type=c('conditional','temporal'),log_tp=FALSE){
  DE_res<-slot(object,'DE_results')
  count_matrix<-exp_matrix(object,'norm')
  for(DE in DE_type){
    if(!(DE %in% names(DE_res))){
      message(paste0(DE,' not found, removing from genes of interest analysis'))
      DE_type<-DE_type[!DE_type==DE]
    }
  }
  list_interest<-list()
  for (DE in DE_type){
    if (DE %in% names(DE_res)){
      DE_list<-DE_res[[DE]]
    }
    for (exp in names(DE_list)){
      sig_df<-DE_list[[exp]][['DE_raw_data']]
      sig_df<-sig_df[sig_df$gene_id %in% genes_of_interest,]
      sig_df<-sig_df[,c('gene_id','log2FoldChange','pvalue','padj')]
      list_interest[[exp]]<-sig_df
    }
  }

  for (gene in genes_of_interest){
    gene_df<-data.frame(NULL)
    for (exp in names(list_interest)){
      sub_gene_df<-list_interest[[exp]]
      if (gene %in% sub_gene_df$gene_id){
        sub_gene_df<-sub_gene_df[sub_gene_df$gene_id==gene,]
        sub_gene_df$exp_name<-exp
        sub_gene_df<-sub_gene_df[,c('exp_name','log2FoldChange','padj','pvalue')]
      }else{
        sub_gene_df<-data.frame(exp_name=exp,log2FoldChange=NA,padj=NA,pvalue=NA)
      }
      if (nrow(gene_df)==0){
        gene_df<-sub_gene_df
      }else{
        gene_df<-rbind(gene_df,sub_gene_df)
      }
    }
    #Only write the dataframe to csv if at least one experiment has values for the gene
    if(all(is.na(gene_df[,'padj'])) == FALSE){
      write.csv(gene_df,paste0(save_location,gene,'.csv'),row.names=FALSE)
    }
    #Plot gene trajectory if gene is in time object
    if(gene %in% row.names(count_matrix)){
      gene_traj_dta<-calculate_gene_traj_data(object,gene,log_tp)
      gene_traj_plot<-plot_single_gene_traj(gene_traj_dta,slot(object,'group_colors'))
      #May cause LOESS warnings if there are too few datapoints
      ggsave(plot=gene_traj_plot,filename = paste0(save_location,gene,'_trajectory.png'))
    }
  }
  upstream_save_loc<-strsplit(save_location,'/')[[1]][1]
  genes_not_found<-genes_of_interest[!genes_of_interest %in% row.names(count_matrix)]
  if(length(genes_not_found)>0){
    write.csv(genes_not_found,paste0(upstream_save_loc,'/','genes_of_interest_not_found.csv'),row.names=FALSE)
  }
}



# utils functions -------------------

#' @title Identify and remove bad combs
#'
#' @description Hclust cannot handle matrices in which for some pairs of rows and columns,
#' only 1 or fewer shared values are non-NA. This function recurrently
#' identifies the most aggravating column/row, excludes that column/row and checks
#' whether more columns/rows need to be excluded
#'
#'
#' Function taken from github user-slagtermaarten
#' from the ComplexHeatmap issue number 155:
#' https://github.com/jokergoo/ComplexHeatmap/issues/155
#'
#'
#' @param mat Matrix to investigate
#' @param min_shared_fields Minimum number of positions that are not NA in both
#' vectors in order not to flag the vector pair as problematic
#'
#' @return list of problematic combs
#'
#' @examples
#' TS_object<-create_example_object_for_R()
#' TS_object <- normalize_timeSeries_with_deseq2(time_object=TS_object)
#' #Perform conditional differential gene expression analysis
#' TS_object<-conditional_DE_wrapper(TS_object)
#' heat_dta<-create_conditional_heatmap_matrix(TS_object)
#' heat_dta<-prepare_heat_data(heat_dta,log_transform=TRUE)
#' problematic_combs<-identify_problematic_combs(heat_dta[['heat_matrix']])
#'
#' @import GenomicRanges
#'
#' @export
#'
#'
identify_problematic_combs <- function(mat, min_shared_fields = 1) {
  exclude_rows <- NULL
  exclude_cols <- NULL
  stopifnot(is.matrix(mat))

  ## Loop over candidate removals
  for (k in 1:nrow(mat)) {
    candidate_rows <- setdiff(1:nrow(mat), exclude_rows)
    problem_row_combs <- NULL
    for (i in candidate_rows) {
      i_idx <- which(candidate_rows == i)
      for (j in candidate_rows[i_idx:length(candidate_rows)]) {
        if (sum(!is.na(mat[i, ]) & !is.na(mat[j, ])) <= min_shared_fields) {
          problem_row_combs <- rbind(problem_row_combs, c(i, j))
        }
      }
    }
    if (is.null(problem_row_combs)) break
    exclude_rows <- c(exclude_rows,
                      as.integer(names(which.max(table(problem_row_combs)))))
  }

  for (k in 1:ncol(mat)) {
    candidate_cols <- setdiff(1:ncol(mat), exclude_cols)
    problem_col_combs <- NULL
    for (i in candidate_cols) {
      i_idx <- which(candidate_cols == i)
      for (j in candidate_cols[i_idx:length(candidate_cols)]) {
        if (sum(!is.na(mat[, i]) & !is.na(mat[, j])) <= min_shared_fields) {
          problem_col_combs <- rbind(problem_col_combs, c(i, j))
        }
      }
    }
    if (is.null(problem_col_combs)) break
    exclude_cols <- c(exclude_cols,
                      as.integer(names(which.max(table(problem_col_combs)))))
  }

  return(list('row' = exclude_rows, 'column' = exclude_cols))
}
