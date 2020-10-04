# Script for loading occurance data quickly
# last revised 2020July17 tpd

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

#tax_reference<-read.csv("data/AnuranTaxRef_20200708.csv")
head(tax_reference)
for(f in tax_reference$final.taxa){
  sp_data <- Anurans_df %>%  
    filter(final.taxa==f)   # only pull out one taxa (matches i)
write.csv(sp_data, paste0("data/occ_data_used/",
                          f, ".csv"))
}

# checking A. americanus -----
US <- map_data("state")#converts state map to dataframe for mapping with ggplot2
usa <- maps::map("usa") #gets map for subsetting entries. Map of the lower 48 obtained through Package ‘maps’
#version 3.3.0 (Becker et al. 2018) in Program R (R Core Team 2016);
main <- map_data("usa", region=c("main"))#subsets map points to only mainland, no islands

df <- data.frame(main$long, main$lat) #dataframe of longitude and latitude 
usa_l48 <- st_as_sf(df, coords=c("main.long","main.lat"),
                    crs=st_crs("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84")) %>%
  st_combine(.) %>%
  st_cast(., "POLYGON")

Anurans_df %>% filter(final.taxa=='Anaxyrus woodhousii') %>%
  count(infraspecificEpithet)
anax<- Anurans_df %>% 
  filter(final.taxa=='Anaxyrus woodhousii') %>%
  st_as_sf(crs=st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84"),
           coords=c("Longitude", "Latitude"))
ggplot()+geom_sf(data=usa_l48) + 
  geom_sf(data = anax, aes(color=year), alpha=0.2)+
  facet_wrap(~infraspecificEpithet)+
  theme(legend.position = 'top')
