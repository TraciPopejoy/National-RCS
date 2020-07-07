# Script for loading occurance data quickly
# last revised 2020July01 tpd

#load libraries
library(tidyverse); library(lubridate)
library(sf); library(rgeos); library(sp); library(rgdal)
library(cowplot)

#list.files('G:/Shared drives/Anuran Biodiversity Project/SDMs/Occurrence Data/Download & Filtering/')[1:4]
main_dir<-"occ_data" #where the occurrence data resides, one folder, one csv for each taxa
taxa<-gsub('_.*','',list.files(main_dir)) #list the taxa we have occurrence data for


for(q in taxa[1:maxsp]){
  occ_file_name<-list.files(main_dir, pattern=q, full.names = T)
  
  #assigning each csv to the correct G.species object
  assign(paste(substr(q,1,1), #pull genus inital out
               gsub("[A-z ]*\\.","", q), #keep everything after the first space
               sep='.'),
         read.csv(occ_file_name) %>% #read the file in
           dplyr::select(-kingdom,-phylum,-class, #don't need
                         -recordNumber, -individualCount, -catalogNumber) %>% #drop for ease of binding the rows together
           mutate(tax=paste(substr(q,1,1), gsub("[A-z ]*\\.","", q),sep="."))) #add a column to indicate main taxa name
  #check this actually worked correctly
  print(paste(paste(substr(q,1,1), gsub("[A-z ]*\\.","", q), sep='.'), #object name
              unique(get(paste(substr(q,1,1), gsub("[A-z ]*\\.","", q),#species within that object name
                               sep='.'))$species), sep="=="))

}

Anurans_df<-bind_rows(mget(paste(substr(taxa[1:maxsp],1,1), #pull genus inital out
                                 gsub("[A-z ]*\\.","", taxa[1:maxsp]), #keep everything after the period
                                 sep='.')))


nobs_source_table<-Anurans_df %>% group_by(tax, species, source) %>% tally() %>%
  pivot_wider(names_from = source, values_from = n)
write.csv(nobs_source_table, paste0('rcs_results/nobs_source_table_',
                                   format(Sys.Date(), '%Y%m%d'),'.csv'))

# need to ask Chloe &/or Meryl about the following species
double_sp<-Anurans_df %>% group_by(tax) %>% 
  summarize(all.spp=paste(unique(species), collapse=", ")) %>%
  group_by(tax, all.spp) %>%
  dplyr::select(tax, all.spp) %>%
  filter(grepl(', ', all.spp))
write.csv(double_sp, paste0('rcs_results/double_sp_',
                                    format(Sys.Date(), '%Y%m%d'),'.csv'))


#temporary fix for double species names
spkey<-Anurans_df %>% group_by(tax) %>% 
  summarize(all.spp=paste(unique(species), collapse=", ")) %>%
  group_by(tax, all.spp) %>%
  dplyr::select(tax, all.spp) %>% 
  mutate(first_sp=gsub(',.*', '', all.spp),
                     second_sp=substring(gsub('.*,', '', all.spp),2)) %>%
  ungroup() %>% dplyr::select(-all.spp) %>%
  gather("df", "species", -tax) %>% dplyr::select(-df)

An_df_temp<-left_join(Anurans_df, spkey)
                     