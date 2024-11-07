#Basic processing libraries
library(dplyr)
library(tidyr)


#Load decouplR libraries
library(decoupleR)
library(OmnipathR)

#Plotting libraries
library(ggplot2)
library(tidyverse)
library(ggseqlogo)
library(ComplexHeatmap)
library(ggrepel)
library(RColorBrewer)


#Load JASPAR libraries
library(JASPAR2024)
library(RSQLite)


# wrapper functions -------------------

#' @title Wrapper for the basic decouplR analysis
#'
#' @description A wrapper which can run both versions of the decouplR pipeline
#' decouplR can be run using counts or degs. Where in the case of counts the entire
#' count matrix of the TimeSeries object is loaded and transcription factors are
#' found based on the amount of variance between the samples.
#'
#' In the DEG format, each differentiatl gene expression analysis of the object is
#' iterated over and significant transcription factors are found from the statistic used
#' by the differential expression pipeline (wald statistic with DESeq2 and t-value with limma)
#'
#' If a save location is provided the function will save a heatmap for either
#' method (DEG or counts).
#'
#' Regardless of save location, the function returns a list of the identified
#' transcription factors. The number is based on the desired number of top TFs.
#'
#' @param time_object A TimeSeries object
#' @param TF_network The collectTri network
#' @param TF_method Either 'counts' or 'DEG' - indicates the decouplR analysis to perform
#' @param save_location The location in which the heatmaps will be saved
#' @param top_TF An integer for the number of top TFs to identify
#' @param TF_statistic The statistic to use for the decouplR analysis
#' @param minimum_TF_size The minimum size of transcription factor to query for
#'
#' @return A list of identified TFs
#'
#' @export
wrapper_TF_extraction<-function(time_object,TF_network,TF_method,save_location=NULL,top_TF=25,TF_statistic='stat',minimum_TF_size=5){

  found_TFs<-c()

  if(TF_method=='counts'){
    sample_file<-exp_sample_data(time_object)
    count_matrix<-exp_matrix(time_object,'norm')

    TF_count_res<-TF_count_analysis(net=TF_network,counts=count_matrix,target_tfs=top_TF,min_size=minimum_TF_size)
    TF_count_res<-as.data.frame(TF_count_res)

    found_TFs<-c(found_TFs,TF_count_res$source)

    if(is.null(save_location)==FALSE){
      dir.create(save_location)

      samp_file<-exp_sample_data(time_object)
      samples<-row.names(TF_count_res)
      new_names<-paste0(samp_file[samples,]$sample,'_',samp_file[samples,]$group,'_',samp_file[samples,]$timepoint)
      row.names(TF_count_res)=new_names

      my_heat<-plot_standard_TF_heatmap(TF_count_res)
      save_pheatmap_pdf(my_heat,filename=paste0(save_location,'top25_tf_count_heat.pdf'))
    }

  }else if(TF_method=='DEG'){
    DE_types<-names(time_object@DE_results)

    DE_lst<-create_DE_list_from_time_object(time_object)
    #This runs the decouplR analysis for each DEG in the object and extracts the top n TFs for each
    TF_df<-TiSA_deg_TF_analysis(DE_lst,TF_network,top_TF,TF_statistic=TF_statistic,minimum_TF_size=minimum_TF_size)

    found_TFs<-unique(TF_df$source)
    if(is.null(save_location)==FALSE){
      dir.create(save_location)

      #We perform another decouplR analysis with a list of TFs this time to get scores and
      #pvalues for each of the found TFs in the topn indiv analysis
      TF_full_df<-TiSA_deg_TF_analysis(DE_lst=DE_lst,TF_network=TF_network,target_TFs=found_TFs,
                                       TF_statistic=TF_statistic,minimum_TF_size=minimum_TF_size)

      #Create heatmap
      create_deg_TF_heatmap(DE_lst,TF_full_df,custom_group_names,custom_colors,
                            save_name=paste0(save_location,'Summary_heat_map.pdf'))
    }
  }else{
    message("Please use 'counts' or 'DE' as the TF method")
    return(NULL)
  }
  return(unique(found_TFs))
}


#' @title A small wrapper to perform the TF analysis for time series objects
#'
#' @description Iterates over a DE list to run the TF decouplR analysis
#'
#' @param DE_lst A list containing the DE files in the TimeSeries Object
#' @param TF_network The collectTri network
#' @param target_TFs Either top n or a list of TFs
#' @param TF_statistic The statistic to use in the decouplR analysis
#' @param minimum_TF_size The minimum size of TFs to allow to be queried.
#'
#' @return The results from the TF analysis
#'
#' @export
TiSA_deg_TF_analysis<-function(DE_lst,TF_network,target_TFs,TF_statistic='stat',minimum_TF_size=5){
  #Create a df with the found TFs to create the heatmap
  TF_full_df<-data.frame(NULL)
  for(exp in names(DE_lst)){
    TF_deg_res<-TF_deg_analysis(net=TF_network,deg=DE_lst[[exp]],target_tfs=target_TFs,stat_type=TF_statistic,min_size=minimum_TF_size)
    TF_deg_res$group<-rep(exp,nrow(TF_deg_res))
    if(nrow(TF_full_df)==0){
      TF_full_df<-TF_deg_res
    }else{
      TF_full_df<-rbind(TF_full_df,TF_deg_res)
    }
  }
  return(TF_full_df)
}

# Processing functions -------------------


#' @title Perform the TF analysis using counts
#'
#' @description Utilizes a count matrix and determines which transcription
#' factors are activated/deactivated in each sample of the count file
#'
#' @param net The CollectTri network
#' @param counts the count data
#' @param target_tfs Can be a vector of TFs or a number for topn
#' @param min_size Minimum number of genes contained in an accepted TF
#'
#' @return A file with the requested number of top TFs from the analysis
#'
#'
#' @export
TF_count_analysis<-function(net,counts,target_tfs,min_size=5,scale=TRUE){
  sample_acts <- run_ulm(mat=counts, net=net, .source='source', .target='target',
                         .mor='mor', minsize = min_size)

  # Transform to wide matrix
  sample_acts_mat <- sample_acts %>%
    pivot_wider(id_cols = 'condition', names_from = 'source',
                values_from = 'score') %>%
    column_to_rownames('condition') %>%
    as.matrix()

  if(is.numeric(target_tfs)==TRUE){
    tfs <- sample_acts %>%
      group_by(source) %>%
      summarise(std = sd(score)) %>%
      arrange(-abs(std)) %>%
      head(target_tfs) %>%
      pull(source)
  }else{
    tfs<-target_tfs
  }

  #Filter per top TFs
  sample_acts_mat <- sample_acts_mat[,tfs]

  # Scale per sample - changes score
  if(scale==TRUE){
    sample_acts_mat <- scale(sample_acts_mat)
  }


  return(sample_acts_mat)
}

#' @title Perform the TF analysis using DEGs from DESeq2
#'
#' @description Function which uses the WALD statistic from DESeq2 to calculate
#' the transcription factor score. Ultimately stating if TFs are activated/supressed
#' in the conditions of the differential gene expression. This analysis also enables the
#' viewing of genes within TFs and their estimated role, as in do they serve to activate
#' a TF or deactivate it.
#'
#' @param net The CollectTri network
#' @param deg the differential gene expression data
#' @param target_tfs Can be a vector of TFs or a number for topn
#' @param stat_type The statistic contained in the data. Usually stat (wald) or t-test
#' @param min_size Minimum number of genes contained in an accepted TF
#'
#' @return The resulting TF data
#'
#'
#' @export
TF_deg_analysis<-function(net,deg,target_tfs,stat_type='stat',min_size=5,do_rank=TRUE){

  #Use the DEG statistic - shows more the conditional aspect rather than the per sample
  contrast_acts <- run_ulm(mat=deg[, stat_type, drop=FALSE], net=net, .source='source', .target='target',
                           .mor='mor', minsize = min_size)

  if(do_rank==TRUE){
    f_contrast_acts<-perform_ranking_for_tf_deg(dta_to_rank=contrast_acts)
  }else{
    f_contrast_acts<-contrast_acts
  }
  #This arranges things in order of rank
  #Since we created the ranks separately it will
  #select an even number of 'top' ranks from positive and negative
  #scores
  if(is.numeric(target_tfs)==TRUE){
    tfs <- f_contrast_acts %>%
      arrange(rnk) %>%
      head(target_tfs) %>%
      pull(source)
  }else{
    tfs<-target_tfs
  }


  #Simply filters for the selected top-ranks from earlier
  f_contrast_acts <- f_contrast_acts %>%
    filter(source %in% tfs)
  return(f_contrast_acts)
}


#' @title Ranks the deg data
#'
#' @description This function ranks the deg data using the score as calculated
#' by the decouplr pipeline. The rank is calculated separetely for degs
#' with negative and positive scores. Rank implies a numerical rank such as
#' 1,2,3 etc... where 1 is the highest score (or lowest if negative scores)
#'
#' @param dta_to_rank The data as produced by \code{TF_deg_analysis}
#'
#' @return The ranked deg data
#'
#' @export
perform_ranking_for_tf_deg<-function(dta_to_rank){
  # Filter top TFs in both signs
  f_contrast_acts <- dta_to_rank %>%
    mutate(rnk = NA)
  msk <- f_contrast_acts$score > 0
  #Rank sorts them and assigns them their 'rank' from largest to smallest
  #Here we rank scores which are positive and negative separately.
  f_contrast_acts[msk, 'rnk'] <- rank(-f_contrast_acts[msk, 'score'])
  f_contrast_acts[!msk, 'rnk'] <- rank(-abs(f_contrast_acts[!msk, 'score']))

  return(f_contrast_acts)

}


#' @title Format the frequency matrix
#'
#' @description Basic function which creates a wide format of the frequency matrix
#' and preserves only the columns representing nucleotides. This is a necessary
#' step before the creation of the sequence matrix to plot seqlogos
#'
#' @param freq_mat Frequency matrix as taken from JASPAR2024 database
#'
#' @return The formatted frequency matrix
#'
#' @export
format_freq_matrix<-function(freq_mat){
  # Convert the data to wide format
  freq_mat <- freq_mat %>%
    pivot_wider(names_from = row, values_from = val, values_fill = 0)

  #Convert to standard R df, limit to necessary columns
  freq_mat<-as.data.frame(freq_mat[,c('A','C','G','T')])
  return(freq_mat)
}

#' @title Calculate nucleotide content
#'
#' @description Calculates the content in percentage of a (or several) nucleotides
#' in a transcription factor. When querying multiple nucleotides the result is given
#' as a single percentage. It is therefore the content of those several nucleotides
#' combined.
#' This is done for the entirety of the JASPAR dataset for the specified species.
#'
#' @param JASPAR_dta The extracted JASPAR dataset
#' @param species_ID The species ID for which the JASPAR dataset will be filtered
#' @param target_nucleotide Vector of nucleotides (or single nucleotide) for which the
#' content will be calculated
#'
#' @return A dataframe containing the relevant JASPAR database columns along with
#' the calculated nucleotide content.
#'
#' @export
calculate_nucleotide_content<-function(JASPAR_dta,species_ID,target_nucleotide=c('A','T')){
  target_species_IDs<-dbReadTable(JASPAR_dta,"MATRIX_SPECIES")$ID[dbReadTable(JASPAR_dta,"MATRIX_SPECIES")$TAX_ID==species_ID]

  main_table <- dbReadTable(db,"MATRIX")
  matrix_data <- dbReadTable(db,"MATRIX_DATA")

  main_table<-main_table[main_table$ID %in% target_species_IDs,]
  matrix_data<-matrix_data[matrix_data$ID %in% target_species_IDs,]

  nucleotide_vect<-c()
  for(TF_ID in main_table$ID){
    frequency_matrix<-matrix_data[matrix_data$ID==TF_ID,]
    nucleotide_counts<-format_freq_matrix(frequency_matrix)

    nucleotide_percentages<-colSums(nucleotide_counts)/(sum(colSums(nucleotide_counts)))*100
    nuc_percent<-sum(nucleotide_percentages[target_nucleotide])

    nucleotide_vect<-c(nucleotide_vect,nuc_percent)
  }

  column_name<-paste0(paste(target_nucleotide,collapse = ''),'_CONTENT')

  main_table[[column_name]]<-nucleotide_vect
  nucleotide_table<-main_table[order(main_table[[column_name]],decreasing = T),]

  return(nucleotide_table)
}

#' @title Create sequence vector
#'
#' @description Function which creates a vector of nucleotide sequences
#' given a formatted frequency matrix.
#'
#' @param nucleotide_counts Formatted frequency matrix
#'
#' @return The sequence vector
#'
#' @export
generate_seq_vect_from_JASPAR_freq_matrix<-function(nucleotide_counts){
  #Calculate number of positions (columns) and possibilities (sum of a row)
  number_of_cols<-nrow(nucleotide_counts)
  number_of_possibilities<-min(rowSums(nucleotide_counts))
  #Create empty vector which will store all DNA sequences
  seq_vector<-c()
  for(number_pos in 1:number_of_possibilities){ #Iterate over possibilities
    #Create empty vector which will store the current iterations DNA sequence
    current_seq<-c()
    for(loc in 1:number_of_cols){
      #Isolate necessary position (called column, represented as rows in the dataframe)
      possible_nucleotides<-nucleotide_counts[loc,]

      #Search for non-zero value nucleotides
      possible_nucleotides<-colnames(possible_nucleotides)[possible_nucleotides>0]

      #randomly select Nucleotide from non-zero options
      chosen_nucleotide<-sample(possible_nucleotides, 1)

      #update order
      current_seq<-c(current_seq,chosen_nucleotide)

      #Update values (-) on the selected nucleotide
      nucleotide_counts[loc,chosen_nucleotide]<-nucleotide_counts[loc,chosen_nucleotide]-1
    }
    #Collapse into the same character
    current_seq<-paste(current_seq,collapse='')
    #Update the main vector
    seq_vector<-c(seq_vector,current_seq)
  }
  return(seq_vector)
}

#' @title Add pvalue to TF data
#'
#' @description Function which performs a decouplr analysis for all TFs. Then function
#' is intended to be run on a large amount of TFs. It goes through all the experiments
#' of a timeseriesobject, appends the found pvalues to a dataframe.
#'
#' The function then filters using the provided pvalue threshold. It returns the filtered
#' table. Along with the combined results from the TF analysis.
#'
#' @param time_object The TimeSeriesAnalysis object
#' @param TF_network The CollectTri network
#' @param nuc_table A table from the JASPAR dataset which also contains nucleotide content
#' This table is produced by \code{calculate_nucleotide_content}
#' @param TF_statisti The statistic contained in the data. Usually stat (wald) or t-test
#' @param minimum_TF_size Minimum number of genes contained in an accepted TF
#' @param filter_pval The pvalue threshold to use on the significance filter
#'
#' @return A list with the nuc_table updated with pvalues and the full TF results.
#'
#' @export
add_pvalue_from_decouplr<-function(time_object,TF_network,nuc_table,TF_statistic='stat',minimum_TF_size=5,filter_pval=0.05){

  TF_df_full<-data.frame(NULL)
  target_TFs<-nuc_table$NAME
  DE_types<-names(time_object@DE_results)
  for(DE in DE_types){
    for(exp in names(time_object@DE_results[[DE]])){
      DE_file<-time_object@DE_results[[DE]][[exp]]$DE_raw_data
      row.names(DE_file)<-DE_file$gene_id
      TF_deg_res<-TF_deg_analysis(net=TF_network,deg=DE_file,target_tfs=target_TFs,stat_type=TF_statistic,min_size=minimum_TF_size)

      for_merge<-TF_deg_res[,c('source','p_value','rnk')]
      colnames(for_merge)<-c('NAME',paste0('p_value_',exp),paste0('rank_',exp))
      nuc_table<-merge(nuc_table,for_merge,by='NAME',all=T)

      TF_deg_res$group<-rep(exp,nrow(TF_deg_res))
      if(nrow(TF_df_full)==0){
        TF_df_full<-TF_deg_res
      }else{
        TF_df_full<-rbind(TF_df_full,TF_deg_res)
      }
    }
  }

  #Identify the required column
  target_col<-colnames(nuc_table)[endsWith(colnames(nuc_table),'CONTENT')]

  #Remove non-JASPAR IDs
  nuc_table <- nuc_table[!is.na(nuc_table$ID), ]
  nuc_table<-nuc_table[order(nuc_table[[target_col]],decreasing = T),]

  #Filter for TFs where any of the pvalue rows is below the threshold
  nuc_table <- nuc_table %>%
    filter(if_any(starts_with("p_value"), ~ . < filter_pval))

  return(list('nucleotide_table'=nuc_table,'TF_full_df'=TF_df_full))
}



#' @title Function which merges TF scores and nucleotide content
#'
#' @description This function merges the information from nucleotide content and
#' transcription factor data. It takes in a dataframe containing nucleotide content
#' and another containing transcription factor scores. It merges this information and
#' adds the nucleotide content score/percentage to the name of the transcription
#' factor. This is mainly done to aid in the interpretation of heatmap plots as they
#' will contain TF names and their associated nucleotide percentage.
#'
#' @param nucleotide_data Dataframe containing the nucleotide content
#' @param TF_data The transcription factor data as produced by decouplR
#'
#' @return None
#'
#' @export
add_nucleotide_percentage_to_TF<-function(nucleotide_data,TF_data){
  #Identify the required column
  target_col<-colnames(nucleotide_data)[endsWith(colnames(nucleotide_data),'CONTENT')]

  nucleotide_data<-nucleotide_data[,c('NAME',target_col)]
  nucleotide_data$NEW_NAME<-paste0(nucleotide_data$NAME,'_',round(nucleotide_data[[target_col]],digits = 1))
  nucleotide_data<-nucleotide_data[,c('NAME','NEW_NAME')]
  colnames(nucleotide_data)=c('source','NEW_NAME')
  nucleotide_data <- nucleotide_data %>%
    distinct(source, .keep_all = TRUE)


  TF_data <- TF_data %>%
    left_join(nucleotide_data %>% select(source, NEW_NAME), by = "source") %>%
    # Replace 'source' with 'NEW_NAME' wherever there's a match
    mutate(source = coalesce(NEW_NAME, source)) %>%
    # Remove the NEW_NAME column since it's no longer needed
    select(-NEW_NAME)

  return(TF_data)
}

# Plotting functions -------------------


#' @title Plots a TF heatmap for counts
#'
#' @description This function uses pheatmap to plot a standard heatmap
#' showing the TF scores (calculated from counts) for each sample
#'
#' @param heat_data The data as produced by \code{TF_count_analysis}
#'
#' @return The pheatmap object
#'
#'
#' @export
plot_standard_TF_heatmap<-function(heat_data){
  # Choose color palette
  palette_length = 100
  my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

  my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
                 seq(0.05, 3, length.out=floor(palette_length/2)))

  # Plot - Measure of variability in the counts
  heat_plot<-pheatmap::pheatmap(as.matrix(heat_data), border_color = NA, color=my_color, breaks = my_breaks)
  trash<-capture.output(dev.off())#Capture output to prevent print
  return(heat_plot)
}


#' @title Save pheatmap
#'
#' @description Simple function to save the grid based plotting of pheatmap. Saves as pdf
#'
#' @param x The pheatmap object
#' @param filename The filename to use
#' @param width The width in cm
#' @param height The height in cm
#'
#' @return None
#'
#'
#' @export
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  trash<-capture.output(dev.off())#Capture output to prevent print
}




#' @title Plot TF DEG barplot results
#'
#' @description Function which creates a barplot of the most activated/deactivated
#' TFs found from the DEG analysis \code{TF_deg_analysis}
#'
#' @param TF_deg_dta data as produced by \code{TF_deg_analysis}
#'
#' @return The ggplot2 object
#'
#'
#' @export
TF_deg_barplot<-function(TF_deg_dta){
  plt<-ggplot(TF_deg_dta, aes(x = reorder(source, score), y = score)) +
    geom_bar(aes(fill = score), stat = "identity") +
    scale_fill_gradient2(low = "darkblue", high = "indianred",
                         mid = "whitesmoke", midpoint = 0) +
    theme_minimal() +
    theme(axis.title = element_text(face = "bold", size = 12),
          axis.text.x =
            element_text(angle = 45, hjust = 1, size =10, face= "bold"),
          axis.text.y = element_text(size =10, face= "bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    xlab("TFs")

  return(plt)
}


#' @title Create a volcano plot from the genes associated to a TF
#'
#' @description Given a specific transcription factor, this function will
#' get the genes involved in that transcription factor. It plots these genes on
#' a volcano plot but utilizes the mor and wald statistic to color them. With
#' a red gene being one that activates the TF and a blue one that deactivates it.
#'
#' A dotted line on the y axis is added to mark the threshold of significance based
#' on the 0.05 standard. The pvalue or adjusted p value is used depending on the content
#' of the provided deg data.
#'
#' @param time_object The TimeSeriesAnalysis object to query for the deg file
#' @param target_TF The target transcription factor
#' @param net The CollectTri network
#' @param target_exp The experiment to take the differential expression dta from
#' @param p_type_used pvalue or padj depending on the statistic that was added to the deg file.
#'
#' @return The ggplot2 object
#'
#'
#' @export
plot_individual_TF<-function(time_object,target_TF,net,DE_type,target_exp,p_type_used='padj',save_location=NULL){

  DE_file<-time_object@DE_results[[DE_type]][[target_exp]]$DE_raw_data
  row.names(DE_file)<-DE_file$gene_id
  if(is.null(DE_file)==TRUE){
    message('No DE file found in the given DE type and target experiment')
    return(NULL)
  }
  df <- net %>%
    filter(source == target_TF) %>%
    arrange(target) %>%
    mutate(ID = target, color = "3") %>%
    column_to_rownames('target')
  inter <- sort(intersect(rownames(DE_file),rownames(df)))
  df <- df[inter, ]
  df[,c('log2FoldChange', 'stat', p_type_used)] <- DE_file[inter, c('log2FoldChange', 'stat', p_type_used)]
  df <- df %>%
    mutate(color = if_else(mor > 0 & stat > 0, 'Activates TF', color)) %>%
    mutate(color = if_else(mor > 0 & stat < 0, 'Deactivates TF', color)) %>%
    mutate(color = if_else(mor < 0 & stat > 0, 'Deactivates TF', color)) %>%
    mutate(color = if_else(mor < 0 & stat < 0, 'Activates TF', color))

  #Here blue means that the sign of multiplying the mor and t-value is negative, meaning that these genes are “deactivating”
  #the TF, and red means that the sign is positive, meaning that these genes are “activating” the TF.

  #Find gene with pvalue (or padj) closest to the 0.05 threshold while remaining smaller
  low_gene_target<-DE_file[DE_file[[p_type_used]]<0.05,]
  low_gene_target<-row.names(low_gene_target)[low_gene_target[[p_type_used]]==max(low_gene_target[[p_type_used]])]
  low_gene_target<-low_gene_target[1]#In case there are multiple of the same value

  #Create a dataframe of necessary columns and convert pvalue to -log10 (as is it's form in the plot)
  find_lowest_sig<-DE_file[,p_type_used]
  find_lowest_sig<- -log10(find_lowest_sig)
  names(find_lowest_sig)=row.names(DE_file)
  #Extract value based on the gene.
  lowest_sig<-unname(find_lowest_sig[low_gene_target])

  highest_value_L2FC<-max(abs(df$log2FoldChange))

  plt<-ggplot(df, aes(x = log2FoldChange, y = -log10(!!sym(p_type_used)), color = color)) +
    geom_point()+
    scale_colour_manual(name='Legend',values = c("red","royalblue3",'grey'),
                        guide = guide_legend(override.aes = list(size = 5))) +
    geom_label_repel(aes(label = ID, size=1),show.legend=FALSE) +
    theme_minimal() +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    #Add a significant dotted line (red)
    geom_hline(yintercept=lowest_sig, linetype="dashed", color = "black")+
    ggtitle(paste0('Transcription factor: ',target_TF,' | Experiment: ',target_exp))+
    xlim(-highest_value_L2FC,highest_value_L2FC)+
    theme(text = element_text(size = 10),
          plot.title = element_text(size = 12),
          legend.text = element_text(size = 10))

  if(is.null(save_location)==FALSE){
    file_name<-paste0(save_location,target_TF,'_',target_exp,'_volcano.png')
    ggsave(plot = plt,filename = file_name,dpi=600,bg='white')
  }else{
    return(plt)
  }

}


#' @title Create the DEG heatmap
#'
#' @description Creates a heatmap which summarizes the various transcription
#' factors identified in the different DE files used prior to this function
#'
#' The function requries a deg dataframe with at least two groups as well as a
#' named vector giving custom group names (with names being the original and
#' values being the new names).
#'
#' Colors must also be provided via named vector, where the names are the name
#' of the groups (use new names) and values are the colors.
#'
#' save_name should contain the full path, name of the file, and .pdf extension.
#'
#' The heatmap is directly saved to the computer, not returned.
#'
#' Note the following on the 'mor' statistic
#' The Mode of Regulation (MoR) indicates whether the inferred activity of a
#' transcription factor (TF) is likely to be activating or repressing.
#' It helps to understand how a transcription factor influences its target genes.
#'
#' @param heat_data The heatmap data
#' @param custom_group_name Named vector where names are old names and values new names
#' @param custom_colors Named vector where names are group names and values are colors
#' @param save_name Full path, file name, and .pdf extension
#'
#' @return None
#'
#' @export
create_deg_TF_heatmap<-function(DE_list,heat_data,custom_group_names,custom_colors,
                                save_name=NULL){
  #Set up custom group names for heatmaps
  custom_group_names<-names(DE_list)
  names(custom_group_names)<-names(DE_list)

  #set up custom colors for the groups - use new groups
  custom_colors<-get_qualitative_colors(length(custom_group_names))
  names(custom_colors)=custom_group_names

  custom_width<-3*length(custom_group_names)
  custom_height<-0.0375*nrow(heat_data)

  data<-heat_data[,c('source','score','group','p_value')]

  #Format the data
  formatted_data<-NULL
  pval_dta<-NULL
  for(i in unique(data$group)){
    temp_data<-data[data$group==i,]

    temp_p_data<-as.data.frame(temp_data[,c('source','p_value')])
    temp_data<-as.data.frame(temp_data[,c('source','score')])

    colnames(temp_p_data)<-c('TFs',i)
    colnames(temp_data)<-c('TFs',i)
    if(is.null(formatted_data)==TRUE){
      formatted_data<-temp_data
      pval_dta<-temp_p_data
    }else{
      formatted_data<-merge(formatted_data,temp_data,by='TFs',all=TRUE)
      pval_dta<-merge(pval_dta,temp_p_data,by='TFs',all=TRUE)
    }
  }
  row.names(formatted_data)=formatted_data$TFs
  formatted_data<-formatted_data[,c(2:ncol(formatted_data))]

  #Replace names
  if(is.null(custom_group_names)==FALSE){
    colnames(formatted_data)=unname(custom_group_names[colnames(formatted_data)])
  }
  #Create coloring map
  color_map=unname(custom_colors[colnames(formatted_data)])

  #Set up the column splits
  split = rep(1:length(colnames(formatted_data)), each = )
  ha = HeatmapAnnotation(
    foo = anno_block(gp = gpar(fill = color_map), labels = colnames(formatted_data))
  )

  #Handle pvalue matrix
  row.names(pval_dta)=pval_dta$TFs
  pval_dta<-pval_dta[,c(2:ncol(pval_dta))]

  # Create a binary matrix to hold '*' for significant cells
  signif_matrix <- matrix("", nrow = nrow(formatted_data), ncol = ncol(formatted_data))
  row.names(signif_matrix) <- row.names(formatted_data)
  colnames(signif_matrix) <- colnames(formatted_data)
  # Populate with '*' where p-values are below 0.05
  signif_matrix[pval_dta < 0.05] <- "*"

  pdf(save_name,height=custom_height,width=custom_width)
  #Plot the heatmap
  draw(Heatmap(as.matrix(formatted_data), name = "TF_score",column_title = NULL,
               column_split = split,top_annotation = ha, show_column_names = FALSE,
               cluster_column_slices = F,
               cell_fun = function(j, i, x, y, custom_width, custom_height, fill) {
                 grid.text(signif_matrix[i, j], x = x, y = y, gp = gpar(col = "black", fontsize = 10))
               })
  )
  trash<-capture.output(dev.off())#Capture output to prevent print
}

#' @title Function that creates and saves seq logos
#'
#' @description Function which retrieves the frequency matrix of a transcription
#' factor, formats it, creates a seq vector (vector of sequences) from it, then finally
#' creates a seq logo plot for the transcription factor.
#' Note that these should be identical to the ones found in the JASPAR database.
#'
#' If a save_location is provided the seqlogo is saved, otherwise it is
#' returned.
#'
#' @param JASPAR_dta The extracted JASPAR dataset
#' @param target_TF_ID The identifier for a TF
#' @param save_location The save location for the seqlogo
#'
#' @return None
#'
#' @export
plot_seqlogo<-function(JASPAR_dta,target_TF_ID,save_location=NULL){

  matrix_data <- dbReadTable(db,"MATRIX_DATA")
  matrix_data<-matrix_data[matrix_data$ID == target_TF_ID,]
  matrix_data<-format_freq_matrix(matrix_data)

  main_table <- dbReadTable(db,"MATRIX")
  TF_name<-main_table$NAME[main_table$ID==target_TF_ID][1]

  freq_vect<-generate_seq_vect_from_JASPAR_freq_matrix(matrix_data)
  custom_width=nchar(freq_vect[1])+1
  seqlogo<-ggseqlogo(freq_vect)
  seqlogo<-seqlogo+ggtitle(TF_name)
  if(is.null(save_location)==FALSE){
    save_name<-paste0(save_location,TF_name,'_seqlogo.png')
    ggsave(save_name,dpi=600,bg='white',height=5,width=custom_width)
  }else{
    return(seqlogo)
  }

}


# Utility functions -------------------


#' @title Get qualitative colors
#'
#' @description Simple function to create a vector of colors based on a given number
#' of required color and a defined palet.
#'
#' @param n Number of colors to retrieve
#' @param set Palett to use
#'
#' @return The color vector
#'
#' @export
get_qualitative_colors <- function(n,set='Set3') {
  # Choose a qualitative color palette
  palette_name <- set  # Other options: "Set1", "Dark2", etc.

  # Get maximum number of colors in the chosen palette
  max_colors <- brewer.pal.info[palette_name, "maxcolors"]

  # Generate colors based on n
  if (n <= max_colors) {
    colors <- brewer.pal(n, palette_name)
  } else {
    # If n is greater than the palette's max, interpolate more colors
    colors <- colorRampPalette(brewer.pal(max_colors, palette_name))(n)
  }

  return(colors)
}


#' @title Create list of DE files
#'
#' @description Function which creates a list of differential expression files
#' where the name is the experiment name and the value is a DE file
#'
#' @param time_object A TimeSeries Object
#'
#' @return The resulting list of DE files
#'
#' @export
create_DE_list_from_time_object<-function(time_object){
  DE_lst<-list()
  DE_types<-names(time_object@DE_results)
  for(DE in DE_types){
    for(exp in names(time_object@DE_results[[DE]])){
      DE_file<-time_object@DE_results[[DE]][[exp]]$DE_raw_data
      row.names(DE_file)<-DE_file$gene_id
      DE_lst[[exp]]<-DE_file
    }
  }
  return(DE_lst)
}




