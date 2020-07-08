# Script for loading occurance data quickly
# last revised 2020July01 tpd

#load libraries
library(tidyverse); library(lubridate)
library(sf); library(rgeos); library(sp); library(rgdal)
library(cowplot)

#list.files('G:/Shared drives/Anuran Biodiversity Project/SDMs/Occurrence Data/Download & Filtering/')[1:4]
main_dir<-"data/occ_data_raw" #where the occurrence data resides, one folder, one csv for each taxa
taxa<-gsub('_.*','',list.files(main_dir)) #list the taxa we have occurrence data for
tax_reference<-read.csv("data/AnuranTaxRef_20200708.csv")

for(q in taxa){
  occ_file_name<-list.files(main_dir, pattern=q, full.names = T)
  
  #assigning each csv to the correct G.species object
  assign(paste(substr(q,1,1), #pull genus inital out
               gsub("[A-z ]*\\.","", q), #keep everything after the first space
               sep='.'),
         read.csv(occ_file_name) %>% #read the file in
           dplyr::select(-recordNumber, -individualCount, -catalogNumber) %>% #drop for ease of binding the rows together
           mutate(tax=paste(substr(q,1,1), gsub("[A-z ]*\\.","", q),sep="."))) #add a column to indicate main taxa name
  #check this actually worked correctly
  print(paste(paste(substr(q,1,1), gsub("[A-z ]*\\.","", q), sep='.'), #object name
              unique(get(paste(substr(q,1,1), gsub("[A-z ]*\\.","", q),#species within that object name
                               sep='.'))$species), sep="=="))

}

Anurans_df<-bind_rows(mget(paste(substr(taxa,1,1), #pull genus inital out
                                 gsub("[A-z ]*\\.","", taxa), #keep everything after the period
                                 sep='.'))) %>%
  left_join(tax_reference)
rm(list=paste(substr(taxa,1,1), #pull genus inital out
              gsub("[A-z ]*\\.","", taxa), #keep everything after the period
              sep='.'))

nobs_source_table<-Anurans_df %>% group_by(final.taxa, source) %>% tally() %>%
  pivot_wider(names_from = source, values_from = n)
write.csv(nobs_source_table, paste0('rcs_results/nobs_source_table_',
                                   format(Sys.Date(), '%Y%m%d'),'.csv'))


# checking the taxa I use are all valid
library(ritis)
testTax<-function(taxa){
  vec<-NULL
  for(w in 1:length(taxa)){
    vec[w]<-grep("invalid",
                 itis_search(q=paste0('nameWOInd:',gsub(' ','\\ ', taxa[w], fixed=T)))$usage,
                 invert=T)
  }
  return(which(vec==0))
}
#synonym_names(itis_search(q="nameWOInd:Dryophytes\\ andersonii")$tsn)