## code to prepare `DATASET` dataset goes here

# Load raw data from .csv file
list_counts<-list()
for(counts in list.files('data-raw/PBMC_example/raw_counts_TS')){
  list_counts[[counts]]<-read.delim(paste0('data-raw/PBMC_example/raw_counts_TS/',counts),header = F)
}
sample_data<-read.csv("data-raw/PBMC_example/sample_file.csv")

AID_TS_data<-list('sample_dta'=sample_data,'counts'=list_counts)
# Save the cleaned data in the required R package location
usethis::use_data(AID_TS_data, overwrite = TRUE)
