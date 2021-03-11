#TPD combining herpmapper & gbif data for v3
#Last altered 08July2020
library(tidyverse);library(lubridate)

#read in each dataset and merge by species 
main_dir<-'G:/Shared drives/Anuran Biodiversity Project/SDMs/Occurrence Data/Download & Filtering/'
AnuranSp_all<-grep(list.files(main_dir), pattern=".doc", invert=T, value=T)
list.files("G:/My Drive/Questions about Anuran Occurrence.gsheet")

for(k in AnuranSp_all[-c(26,53,93,98)]){
  gbif_folder<-paste(k, list.files(paste0(main_dir, k), pattern="GBIF"), sep='/')
  gbif_file<-list.files(paste0(main_dir, gbif_folder), #makes the directory for each sp
             pattern="v6.csv",#only call the last csv
             full.names = T) #print the entire name so I avoid paste monstrocities
  if(length(gbif_file)==0){GS_GBIF<-data.frame("Latitude"=numeric(),"Longitude"=numeric(),"year"=numeric(),"source"=character(), 
                                           "species"=character(),"infraspecificEpithet"=character(),
                                           "recordNumber"=character(), "individualCount"=numeric(), "catalogNumber"=character())}
  else{
  GS_GBIF <- read.csv(gbif_file) %>% 
    rename("Latitude"="decimalLatitude","Longitude"="decimalLongitude") %>% #match HM dataframe
    mutate(source="GBIF",  #add the source
           infraspecificEpithet=as.character(infraspecificEpithet)) #forcing this to be character for easier merge below
  }
  
  hm_folder<-paste(k, list.files(paste0(main_dir, k), pattern="HerpMapper"), sep='/')
  hm_file<-list.files(paste0(main_dir, hm_folder), #makes the directory for each sp
                        pattern="v3.csv",#only call v3 csv (dates still intact)
                        full.names = T)
  
  if(length(hm_file)==0){GS_HM<-data.frame("Latitude"=numeric(),"Longitude"=numeric(),"year"=numeric(),"source"=character(), 
                                                 "species"=character(),"infraspecificEpithet"=character())}
  else{
    GS_HM <- read.csv(hm_file) %>%
      mutate(source="HM", #add the source
           #new problem: what do to if year doesn't make sense (i.e. 0000-00-00)
           #for now, putting as 2015-01-01 or Sys.Date.... ###### need to ask Meryl & Chloe
           #problem: some records have yyyy-mm-00 as the date (UUID: 4616f3d4-ced4-102d-80d2-00142272a11d)
           #this means as.Date can't parse it as a Date
           #solution: add one to the day to make it valid, this allows it to parse correctly
           DateF=case_when(substr(Date,1,4)=="0000"~as.character(Sys.Date()),
                           substr(Date,9,10)=="00" ~ paste0(substr(Date,1,8),
                                                            as.numeric(substr(Date,9,10))+1),
                           T~Date),
           Date=as.Date(DateF, format="%Y-%m-%d"), #ensure Date column is in Date format
           year=year(Date)) %>%  #pull year out from Date column
      dplyr::select(-DateF) %>% #remove column used to create correct dates (can't have 00)
      rowwise() %>%
      mutate(species=paste(scan(text=Taxon, what=" ", n=2, quiet=T), collapse=" "), #paste the species
             infraspecificEpithet=scan(text=Taxon, what=" ", quiet=T)[3])%>% #paste the subspecies
      ungroup() %>%
      filter(!duplicated(.[c('Latitude', 'Longitude', 'year')])) #remove the duplicates (sensu CEM)
  }
  #combine GBIF and HM datasets by year, lat, long, and source columns
  merged <- full_join(GS_GBIF, GS_HM,
                    by=c("Latitude","Longitude","year","source", 
                         "species","infraspecificEpithet")) %>%
    filter(!duplicated(.[c('Latitude', 'Longitude','year')]),
           !is.na('Latitude'))
  
  output_file<-paste0("occ_data/", gsub(' ','.',k),
                     "_merged_", format(Sys.Date(), "%y%m%d"), "_v3.csv")
  
  write.csv(merged, file = output_file)
  print(paste(k,unique(merged$species), sep='='))
}



#### when you get an:
# Error in file(file, "rt") : invalid 'description' argument
# it means one of the file names isn't working. 
# typically occurs for herp mapper sheets
# species that are excluded from this for loop:
AnuranSp_all[c(26,53,93,98)]
#figure out where the error is occuring:
k
unique(GS_GBIF$species); unique(GS_HM$species)
which(AnuranSp_all==k)


# loading Rana virgatipes by hand - recent synonymized and has two GBIF folders
RV_gbif_file<-list.files(paste0(main_dir, AnuranSp_all[98], "/", AnuranSp_all[98], " GBIF"), #makes the directory for each sp
                        pattern="v6.csv",#only call the last csv
                        full.names = T) #print the entire name so I avoid paste monstrocities
RV_GBIF <- read.csv(RV_gbif_file) %>% 
  rename("Latitude"="decimalLatitude","Longitude"="decimalLongitude") %>% #match HM dataframe
  mutate(source="GBIF",  #add the source
         infraspecificEpithet=as.character(infraspecificEpithet)) 
LV_gbif_file<-list.files(paste0(main_dir, AnuranSp_all[98], "/Lithobates virgatipes GBIF"), #makes the directory for each sp
                         pattern="v6.csv",#only call the last csv
                         full.names = T) #print the entire name so I avoid paste monstrocities
LV_GBIF <- read.csv(LV_gbif_file) %>% 
  rename("Latitude"="decimalLatitude","Longitude"="decimalLongitude") %>% #match HM dataframe
  mutate(source="GBIF",  #add the source
         infraspecificEpithet=as.character(infraspecificEpithet)) 

RV_hm_file<-list.files(paste0(main_dir, AnuranSp_all[98], "/", AnuranSp_all[98], " HerpMapper"), #makes the directory for each sp
                      pattern="v3.csv",#only call v3 csv (dates still intact)
                      full.names = T)
RV_HM <- read.csv(RV_hm_file) %>%
      mutate(source="HM", #add the source
             #new problem: what do to if year doesn't make sense (i.e. 0000-00-00)
             #for now, putting as 2015-01-01 or Sys.Date.... ###### need to ask Meryl & Chloe
             #problem: some records have yyyy-mm-00 as the date (UUID: 4616f3d4-ced4-102d-80d2-00142272a11d)
             #this means as.Date can't parse it as a Date
             #solution: add one to the day to make it valid, this allows it to parse correctly
             DateF=case_when(substr(Date,1,4)=="0000"~as.character(Sys.Date()),
                             substr(Date,9,10)=="00" ~ paste0(substr(Date,1,8),
                                                              as.numeric(substr(Date,9,10))+1),
                             T~Date),
             Date=as.Date(DateF, format="%Y-%m-%d"), #ensure Date column is in Date format
             year=year(Date)) %>%  #pull year out from Date column
      dplyr::select(-DateF) %>% #remove column used to create correct dates (can't have 00)
      rowwise() %>%
      mutate(species=paste(scan(text=Taxon, what=" ", n=2, quiet=T), collapse=" "), #paste the species
             infraspecificEpithet=scan(text=Taxon, what=" ", quiet=T)[3])%>% #paste the subspecies
      ungroup() %>%
      filter(!duplicated(.[c('Latitude', 'Longitude', 'year')])) #remove the duplicates (sensu CEM)
  #combine GBIF and HM datasets by year, lat, long, and source columns
merged <- bind_rows(RV_GBIF, LV_GBIF) %>% full_join(RV_HM,
                      by=c("Latitude","Longitude","year","source", 
                           "species","infraspecificEpithet")) %>%
    filter(!duplicated(.[c('Latitude', 'Longitude','year')]),
           !is.na('Latitude')) 
write.csv(merged, file = paste0("occ_data/", gsub(' ','.',AnuranSp_all[98]),
                                "_merged_", format(Sys.Date(), "%y%m%d"), "_v3.csv"))
