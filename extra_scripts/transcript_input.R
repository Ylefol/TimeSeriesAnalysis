# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#
# BiocManager::install("tximport")
# BiocManager::install("tximportData")
# BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
# BiocManager::install("rhdf5")


library(tximportData)
library(DESeq2)
library(readr)
library(tximport)
library(SummarizedExperiment)
library(TimeSeriesAnalysis)

dir <- system.file("extdata", package = "tximportData")
list.files(dir)
samples <- read.table(file.path(dir, "samples.txt"), header = TRUE)
gene_conv <- read_csv(file.path(dir, "tx2gene.gencode.v27.csv"))
files <- file.path(dir, "kallisto", samples$run, "abundance.tsv.gz")
names(files) <- samples$sample




mock_sample_data<-data.frame(sample=samples$sample,
                             group=c('control','control','control','control','treat','treat'),
                             replicate=c('C_1','C_2','C_1','C_2','T_1','T_1'),
                             timepoint=c(1,1,5,5,1,5))
TS_object <- new('TimeSeries_Object',
                 group_names=c('treat','control'),group_colors=c("#e31a1c","#1f78b4"),DE_method='DESeq2',
                 DE_p_filter='padj',DE_p_thresh=0.05,DE_l2fc_thresh=1,
                 PART_l2fc_thresh=1,sem_sim_org='org.Hs.eg.db',Gpro_org='hsapiens')



add_exp_data_kallisto<-function(time_object,sample_data,kallisto_files,tx2gene,matrix_name='raw'){
  txi.kallisto <- tximport(kallisto_files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)
  condition<-factor(sample_data$group,levels=rev(time_object@group_names))

  col_data <- data.frame(row.names=colnames(kallisto_files),condition)
  dds <- DESeqDataSetFromTximport(txi.kallisto, col_data, ~condition)


  exp_data<-SummarizedExperiment(assays = assays(dds)$counts,colData = sample_data)
  names(assays(exp_data))=matrix_name

  time_object@exp_data<-exp_data
  return(time_object)
}


test_TS<-add_exp_data_kallisto(TS_object,mock_sample_data,files,gene_conv)
test_TS <- normalize_timeSeries_with_deseq2(time_object=test_TS)
#Perform conditional differential gene expression analysis
test_TS<-conditional_DE_wrapper(test_TS)
#Perform temporal differential gene expression analysis
test_TS<-temporal_DE_wrapper(test_TS,do_all_combinations=T)


#Rest of the pipeline will be fine from here
