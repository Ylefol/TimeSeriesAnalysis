# Load the library
library(circlize)
library(TimeSeriesAnalysis)

#Create functions

#Compression function
compress <- function(values, max_val, max_compression) {
  midpoint <- max_val / 2
  compression_factor <- max_compression * abs(values - midpoint) / (max_val / 2)
  compressed_values <- values - compression_factor * (values - midpoint)
  return(compressed_values)
}

#Load data from TiSA object
load_data_from_TiSA<-function(time_object,gtf_data,chrom_order=c(paste0("chr", 1:22), "chrX", "chrY"),l2fc_colors=c('blue','red')){
  dta_lst<-list()
  for(exp in names(time_object@DE_results$conditional)){
    temp_df<-time_object@DE_results$conditional[[exp]]$DE_sig_data
    temp_df<-temp_df[temp_df$gene_id %in% target_genes,c('gene_id','log2FoldChange','padj')]
    #Get gtf_data
    gtf_data<-gtf_data[,c(1:4)]
    temp_df<-merge(gtf_data,temp_df,by='gene_id')
    temp_df<-temp_df[,c(2:6)]
    temp_df<-temp_df[temp_df$seqid %in% chrom_order,]
    temp_df$seqid<-factor(temp_df$seqid,levels=chrom_order)
    # Order the dataframe
    temp_df <- temp_df[order(temp_df$seqid, temp_df$start), ]
    
    #Set colors
    temp_df$color<-rep(c('black'),nrow(temp_df))
    temp_df$color[temp_df$log2FoldChange<0]<-l2fc_colors[1]
    temp_df$color[temp_df$log2FoldChange>0]<-l2fc_colors[2]
    
    dta_lst[[exp]]<-temp_df
  }
  return(dta_lst)
}

#Compute size for a shared scale
compute_size_same_scale<-function(data_list,target_column='padj'){
  # Get the global min and max of the 'padj' column across all dataframes
  global_min <- min(sapply(data_list, function(df) min(df[[target_column]], na.rm = TRUE)))
  global_max <- max(sapply(data_list, function(df) max(df[[target_column]], na.rm = TRUE)))
  
  # Scale the target_column column for each dataframe
  data_list <- lapply(data_list, function(df) {
    df$size <- 0.3 + (df[[target_column]] - global_min) * (1 - 0.3) / (global_max - global_min)
    return(df)
  })
  return(data_list)
}


#Compress the data to avoid gene points exiting the boundaries
compress_gene_locations<-function(data_list,cytoband_dta,max_compression=0.05){
  for(chrom in unique(data_list[[1]]$seqid)){
    max_value=tail(cytoband_dta$V3[cytoband_dta$V1==chrom],n=1)
    for(exp in names(data_list)){
      data_list[[exp]][data_list[[exp]]$seqid==chrom,]$start<-compress(data_list[[exp]][data_list[[exp]]$seqid==chrom,]$start, max_value, max_compression)
      data_list[[exp]][data_list[[exp]]$seqid==chrom,]$end<-compress(data_list[[exp]][data_list[[exp]]$seqid==chrom,]$end, max_value, max_compression)
    }
  }
  return(data_list)
}

#Create files for genmoic lines
create_genomic_lines_data<-function(gtf_data,cytoband_dta=NULL,max_compression=0.05){
  
  if(is.null(cytoband_dta)==FALSE){
    #put gtf in a list for data format compatibility
    gtf_list<-list(dummy=gtf_data)
    gtf_list<-compress_gene_locations(gtf_list,cytoband_dta,max_compression)
    gtf_data<-gtf_list[['dummy']]
  }
  
  #Pathway has 2*10^8 locations
  bp_per_cluster<-2*10^8/length(unique(gtf_data$cluster))
  gtf_data$mock_location<-rep(0,nrow(gtf_data))
  gtf_data<-gtf_data[order(gtf_data$cluster,gtf_data$seqid,gtf_data$start),]
  
  
  for(idx in 1:length(unique(gtf_data$cluster))){
    target_clust<-unique(gtf_data$cluster)[idx]
    if(idx==1){
      start_loc=0
      end_loc=bp_per_cluster
    }else{
      start_loc=bp_per_cluster*(idx-1) + 1
      end_loc=start_loc+bp_per_cluster
    }
    num_iter<-nrow(gtf_data[gtf_data$cluster==target_clust,])
    if(idx==length(unique(gtf_data$cluster))){#Last iteration, add an iter
      # Adding an extra element allows to use the end of the mock chromosome as the end value
      evenly_spaced_values <- seq(start_loc, end_loc, length.out = num_iter+1)
      gtf_data$mock_location[gtf_data$cluster==target_clust]<-evenly_spaced_values[-length(evenly_spaced_values)]
    }else{
      evenly_spaced_values <- seq(start_loc, end_loc, length.out = num_iter)
      gtf_data$mock_location[gtf_data$cluster==target_clust]<-evenly_spaced_values
    }
  }
  
  genomic_lines_file_1<-gtf_data[,c('seqid','start','end')]
  genomic_lines_file_2<-gtf_data[,c('seqid','mock_location')]
  
  #Adjust second  file as it marks the 'end' locations
  mock_end<-c(genomic_lines_file_2$mock_location[-1],2*10^8)
  genomic_lines_file_2$end<-mock_end
  #Modify column names, add necessary data for line location
  colnames(genomic_lines_file_2)=colnames(genomic_lines_file_1)
  genomic_lines_file_2$seqid<-rep('clusters',nrow(genomic_lines_file_2))
  
  #Store in list for returning
  genomic_lines_lst<-list(chrom_data=genomic_lines_file_1,clust_data=genomic_lines_file_2,gtf_mod=gtf_data)
  return(genomic_lines_lst)
}

#Make the circos plot
create_TiSA_circos<-function(data_list,cytoband_data,genomic_lines_list, highlight_chrom=NULL, clust_label_size=0.1,
                             filename='TiSA_circos.png',width_in=8,height_in=8,resolution=600){
  png(filename=filename, units="in", width=width_in, height=height_in, res=resolution)
  circos.initializeWithIdeogram(cytoband)
  
  #Establish limits - balance y axis limits by selecting the largest of the two limits
  toplim = ceiling(max(sapply(data_list, function(x) max(x$log2FoldChange))))
  lowlim = floor(min(sapply(data_list, function(x) min(x$log2FoldChange))))
  dta_lim<-max(toplim,abs(lowlim))
  
  for(exp in rev(names(data_list))){
    #Create data
    exp_dta<-data_list[[exp]]
    mock_bed<-exp_dta[,c('seqid','start','end','log2FoldChange')]
    colnames(mock_bed)<-c('chr','start','end','value1')
    
    circos.genomicTrackPlotRegion(mock_bed, ylim = c(-dta_lim, dta_lim), panel.fun = function(region, value, ...) {
      current_chrom<-get.cell.meta.data('sector.index')
      
      col=exp_dta$color[exp_dta$seqid==current_chrom]
      size=exp_dta$size[exp_dta$seqid==current_chrom]
      
      circos.genomicPoints(region, value, col = col, cex = size, pch = 16)
      cell.xlim = get.cell.meta.data("cell.xlim")
      circos.lines(cell.xlim, c(0, 0), col = "#00000040")
      
    }, track.height = 0.1)
    
  }
  
  heatmap_data<-create_heatmap_for_clusters(genomic_lines_list)

  heatmap_track<-length(names(data_list))+2
  # cell.ylim = get.cell.meta.data("cell.ylim",sector.index = 'clusters',track.index = heatmap_track)
  # 
  # circos.genomicRect(heatmap_data[,c('start','end')], heatmap_data$value, col = heatmap_data$color, 
  #                    border = heatmap_data$colo, lwd = 0.1,
  #                    ybottom=cell.ylim[1], ytop=cell.ylim[2],
  #                    track.index = heatmap_track,sector.index = 'clusters')
  
  circos.genomicLink(genomic_lines_list[['chrom_data']], genomic_lines_list[['clust_data']], col = genomic_lines_list[['gtf_mod']]$cluster_col)
  if(is.null(highlight_chrom)==FALSE){
    highlight.chromosome(highlight_chrom)
  }
  
  #Clear empty cluster tracks
  for(i in 1:(heatmap_track-1)){
    print(i)
    circos.updatePlotRegion(sector.index = "clusters", track.index = i, bg.border = 'white',bg.col='white')
  }
  
  # Hiodding the regions with the above code leaves some faint lines, so we do this to ensure the result is correct.
  #For padding - first value = bottom, second = left, third = right, last = top
  highlight.sector(sector.index = 'clusters',track.index = c(1:(heatmap_track)),col='white',padding=c(0.01,0.01,0.01,0.01))
  
  #Add cluster labeling
  positions <- heatmap_data$label_loc
  labels <- unique(genomic_lines_list[['gtf_mod']]$cluster)
  
  sector_limits <- get.cell.meta.data("xlim", sector.index = "clusters", track.index = 5)
  track_ylims <- get.cell.meta.data("ylim", sector.index = "clusters", track.index = 5)
  
  for (i in seq_along(positions)) {
    circos.text(x = positions[i], y = track_ylims[2] + 0.5, # Place labels slightly outside the track
                labels = labels[i], facing = "bending.inside", niceFacing = TRUE,
                sector.index = "clusters", track.index = 5, cex=clust_label_size, adj = c(0.5, 5.5) # Center the labels
    )
  }
  
  circos.clear()
  dev.off()
}

create_heatmap_for_clusters<-function(genomic_lines_list){
  map_df<-cbind(genomic_lines_list$clust_data,genomic_lines_list$gtf_mod[,c('cluster','cluster_col')])
  colnames(map_df)=c('seqid','start','end','cluster','cluster_col')
  
  heat_df<-data.frame(NULL)
  start_clust<-c()
  end_clust<-c()
  col_vect<-c()
  for(clust in unique(map_df$cluster)){
    start_clust<-c(start_clust,map_df$start[map_df$cluster==clust][1])
    end_clust<-c(end_clust,tail(map_df$end[map_df$cluster==clust],n=1))
    col_vect<-c(col_vect,map_df$cluster_col[map_df$cluster==clust][1])
  }
  values<-sample(x=1:20,size=length(start_clust))
  
  heat_df<-data.frame(start=start_clust,end=end_clust,value=values,color=col_vect)
  heat_df$label_loc<-heat_df$start+((heat_df$end-heat_df$start)/2)
  return(heat_df)
}

#Load the data
load('~/A_Projects/EpiGen/R_Work_Folder/TimeSeriesAnalysis/results_folder/TS_results_PBMC_example_test_script_method/TiSA_object.Rdata')

#Find max
species='hg38'

# Show all PART genes for each conditional results - Best way to show it
target_genes<-row.names(TS_object@PART_results$part_data)

# I'll therefore need the gtf file to create this.
library(rtracklayer)
my_gtf<-readGFF('~/Downloads/gencode.v47.primary_assembly.annotation.gtf')

#Check genes which are not available?
unique(target_genes[!target_genes %in% my_gtf$gene_name])

#Manual corrections to some gene names - change in gtf file to match PART output
my_gtf$gene_name[my_gtf$gene_name=='H3C2']<-'HIST1H3B'
my_gtf$gene_name[my_gtf$gene_name=='H3C8']<-'HIST1H3G'
my_gtf$gene_name[my_gtf$gene_name=='H3C7']<-'HIST1H3F'
my_gtf$gene_name[my_gtf$gene_name=='H1-1']<-'HIST1H1A'
my_gtf$gene_name[my_gtf$gene_name=='SERTAD4BP']<-'LINC00643'
my_gtf$gene_name[my_gtf$gene_name=='GLUD1P2']<-'GLUD1P7'
my_gtf$gene_name[my_gtf$gene_name=='GAS6-DT']<-'GAS6-AS2'
my_gtf$gene_name[my_gtf$gene_name=='ERFE']<-'FAM132B'
my_gtf$gene_name[my_gtf$gene_name=='BMS1P1']<-'BMS1P5'
my_gtf$gene_name[my_gtf$gene_name=='CAVIN2']<-'SDPR'
my_gtf$gene_name[my_gtf$gene_name=='POU2AF3']<-'COLCA2'
my_gtf$gene_name[my_gtf$gene_name=='BMAL2']<-'ARNTL2'
my_gtf$gene_name[my_gtf$gene_name=='ATP13A3-DT']<-'LINC00884'
my_gtf$gene_name[my_gtf$gene_name=='EPOP']<-'C17orf96'
my_gtf$gene_name[my_gtf$gene_name=='JCAD']<-'KIAA1462'
my_gtf$gene_name[my_gtf$gene_name=='DLEU1']<-'DLEU7-AS1'
my_gtf$gene_name[my_gtf$gene_name=='PLAAT1']<-'HRASLS'
my_gtf$gene_name[my_gtf$gene_name=='ADGRF2P']<-'ADGRF2'
my_gtf$gene_name[my_gtf$gene_name=='IRAG1-AS1']<-'MRVI1-AS1'
my_gtf$gene_name[my_gtf$gene_name=='DIPK1B']<-'FAM69B'
my_gtf$gene_name[my_gtf$gene_name=='MTARC2']<-'MARC2'
my_gtf$gene_name[my_gtf$gene_name=='MROCKI']<-'LINC01268'
my_gtf$gene_name[my_gtf$gene_name=='LORICRIN']<-'LOR'
my_gtf$gene_name[my_gtf$gene_name=='TMEM266']<-'C15orf27'
my_gtf$gene_name[my_gtf$gene_name=='JAK3']<-'T'
my_gtf$gene_name[my_gtf$gene_name=='CAVIN3']<-'PRKCDBP'

# Contains the locations we need
my_gtf<-my_gtf[(my_gtf$gene_name %in% target_genes) &(my_gtf$type=='gene'),c('gene_name','seqid','start','end')]
colnames(my_gtf)[colnames(my_gtf)=='gene_name']<-'gene_id'

#Add color to gtf
cmap<-TS_object@PART_results$cluster_map
cmap$gene_id<-row.names(cmap)
my_gtf<-merge(my_gtf,cmap[,c('gene_id','cluster','cluster_col')],by='gene_id')


#Load TiSA data and merge with gtf
dta_lst<-load_data_from_TiSA(time_object = TS_object,gtf_data = my_gtf)

#Compute the size on a shared scale
dta_lst<-compute_size_same_scale(dta_lst)

#Create custom cytoband
cytoband = read.cytoband()$df
cytoband = rbind(cytoband,
                 data.frame(V1 = "clusters", V2 = 1,  V3 = 2e8, V4 = "", V5 = "")
)

#Compress to avoid extreme gene points overlapping with segment boundaries
dta_lst<-compress_gene_locations(data_list = dta_lst,cytoband_dta = cytoband, max_compression = 0.05)

#Create file for genomic lines
genomic_lines_and_gtf_lst<-create_genomic_lines_data(gtf_data = my_gtf,cytoband_dta = cytoband,max_compression = 0.05)

#Create the plot
create_TiSA_circos(data_list = dta_lst,cytoband_data = cytoband,genomic_lines_list = genomic_lines_and_gtf_lst,
                   clust_label_size=0.5, highlight_chrom = 'chr1')








