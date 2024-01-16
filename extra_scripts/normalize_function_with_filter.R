
normalize_timeSeries_with_deseq2_modified_YL <- function(time_object,controlGenes=NULL,min_samples=NULL,min_counts=NULL,use_batch=FALSE,no_condition=FALSE){
  samp_dta_full<-exp_sample_data(time_object)
  count_matrix<-exp_matrix(time_object,'raw')

  #If only the temporal aspect is of interest - case and control must be the same name
  if(no_condition==TRUE){
    condition<-factor(samp_dta_full$timepoint,levels=rev(unique(samp_dta_full$timepoint)))
  }else{
    #Set levels to the inverse of the submitted group names, therefore control is the reference
    condition<-factor(samp_dta_full$group,levels=rev(time_object@group_names))
  }

  my_matrix<-as.matrix(round(count_matrix))

  #Create a coldata frame and instantiate the DESeqDataSet
  if(use_batch==TRUE){
    batch<-factor(samp_dta_full$batch,levels=unique(samp_dta_full$batch))
    col_data <- data.frame(row.names=colnames(my_matrix),condition,batch)
    dds <- DESeqDataSetFromMatrix(countData=as.matrix(my_matrix), colData=col_data, design=~batch+condition)
  }else{
    col_data <- data.frame(row.names=colnames(my_matrix),condition)
    dds <- DESeqDataSetFromMatrix(countData=as.matrix(my_matrix), colData=col_data, design=~condition)
  }

  if(is.null(controlGenes)==TRUE){
    dds = estimateSizeFactors(object=dds)
  }else{
    numerical_control_index=which(rownames(as.matrix(my_matrix)) %in% controlGenes)
    dds = estimateSizeFactors(object=dds,controlGenes=numerical_control_index)#Ability to input custom sequences
  }

  if(is.null(min_counts)==FALSE || is.null(min_samples)==FALSE){
    #Filtering and quick checks
    if(is.null(min_counts)==FALSE & is.null(min_samples)==FALSE){
      keep <- rowSums(counts(dds) >= min_counts) >= min_samples
      dds <- dds[keep,]
    }else if(is.null(min_counts)==FALSE & is.null(min_samples)==TRUE){
      keep <- rowSums(counts(dds)) >= min_counts
      dds <- dds[keep,]
    }else if(is.null(min_counts)==TRUE & is.null(min_samples)==FALSE){
      stop('For filtering, a value must also be given to minimum counts (min_counts)')
    }
    print("Filtering causes the 'raw' matrix to be filtered as well")
    norm_counts <- as.data.frame(counts(dds,normalized=TRUE))
    norm_counts<- na.omit(norm_counts)

    #re_create exp_data raw matrix to reflect the filtered norm
    raw_matrix<-assays(time_object@exp_data)$raw
    raw_matrix<-raw_matrix[row.names(norm_counts),]

    sample_data<-exp_sample_data(time_object)

    exp_data<-SummarizedExperiment(assays = raw_matrix,colData = sample_data)

    names(assays(exp_data))=c('raw')
    time_object@exp_data<-exp_data
  }else{
    norm_counts <- as.data.frame(counts(dds,normalized=TRUE))
    norm_counts<- na.omit(norm_counts)
  }
  assays(time_object@exp_data)$norm<-norm_counts
  time_object@DESeq2_obj<-dds

  return(time_object)
}

