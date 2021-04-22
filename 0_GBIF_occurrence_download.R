library(rgbif); library(tidyverse); library(sf)

# don't pull down occurrences for species I'm not going to use ----
# excluding those with < 20% range in US
sp_range_wiL48 <-read.csv('rcs_results/old_data/sp_range_within_US.csv') %>% 
  dplyr::select(-X)
# removing exotics based on NIAS
NAS_exotic<-c('Osteopilus septentrionalis','Xenopus laevis')
# others not looking into
mis_others<-c('Rhinella marina', #invasive around the world
              'Lithobates fisheri', #taxonomic nightmare
              'Incilius valliceps') #only native to Mexico
sp_exclude <- sp_range_wiL48 %>% filter(per_inside_us <= .200) %>%
  pull(scientific_name) %>% c(NAS_exotic, mis_others)
sp_exclude

# get gbif species keys ----
anuran.gbif.taxa<-NULL
for(u in 1:nrow(sp_range_wiL48)){
  anuran.gbif.taxa1<-name_backbone(name=sp_range_wiL48[u,1], 
                                   rank='species', kingdom='animals')
  anuran.gbif.taxa<-bind_rows(anuran.gbif.taxa1, anuran.gbif.taxa)
}
anuran.taxa.table<-read.csv('./data/AnuranTaxRef_20200708.csv')
anuran.gbif.taxa<-anuran.gbif.taxa %>%
  filter(!(species %in% sp_exclude), #remove species not interested in 
         !duplicated(speciesKey))  #remove duplicated species to not collect twice
  

# build WKT ----
na_poly<-read_sf('/home/tracidubose/RCS_Anuran_input/na_poly_filter/')

library(rgeos)
na_poly_simp<-na_poly %>% 
  as_Spatial() %>% 
  spTransform(CRSobj = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs") %>% 
  gBuffer(width=10) %>% 
  gSimplify(spgeom=., tol=15)
plot(na_poly_simp)

na_poly_wkt<-na_poly_simp %>% spTransform(CRSobj = "+proj=longlat +ellps=WGS84 +datum=WGS84")
library(wicket)
na_land_wkt<-sp_convert(na_poly_wkt)
nchar(na_land_wkt)

saveRDS(na_poly_wkt, file='north_american_poly.rds')
saveRDS(na_land_wkt, file='north_america_WKT.rds')

na_land_wkt<-readRDS('data/north_america_WKT.rds')

# set up a for loop to pull down occ (> 100k occ) ----
# automate pulling from species list
ind_mat<-matrix(c(1,55,56,92), ncol=2, byrow = T)
anuran_occ<-NULL #where all occ kept eventually
gbif_codes<- NULL #keep the download codes

for(j in 1:2){
  gbif_download = occ_download(
    #step 1 - get occ for a species
    pred_in("taxonKey",  anuran.gbif.taxa$speciesKey[ind_mat[j,1]:ind_mat[j,2]]), 
    pred("hasCoordinate", TRUE), #step 2
    pred("hasGeospatialIssue", FALSE), #partial step 3
    pred_not(pred("issue", "IDENTIFIED_DATE_UNLIKELY")), #step 3
    pred_not(pred("issue", "RECORDED_DATE_MISMATCH")), #step 3
    pred_notnull("country"), #step 4
    pred_within(na_land_wkt), #step 5
    pred_lte("eventDate", "2019-12-25"), #last pull
    pred_lt("coordinateUncertaintyInMeters",2000), #removing uncertain pts
    format = "SIMPLE_CSV",
    user=user,pwd=pwd,email=email)
  
  occ_download_wait(gbif_download) 
  occ_download_get(gbif_download[1], path = "./data")
  unzip(paste0("./data/", gbif_download[1],".zip"), exdir="data")
  occ1<-read_tsv(paste0("./data/", gbif_download[1],".csv"),
                 col_types="cccccccccccccccccccicnnnnnnnnTnnnccccccccTcccccTcc")
  
  anuran_occ<-bind_rows(anuran_occ, occ1)
  gbif_codes<-append(gbif_codes, gbif_download[1])
}
write.table(gbif_codes,
            "GBIF_codes_20210330.txt")
nrow(anuran_occ)
91460+62481
View(anuran_occ)

# remove occ with NA
anuran_occ %>% filter(is.na(decimalLatitude)) %>% nrow()

write.csv(anuran_occ, 'anuranOccOrig.csv')

anuran_occ<-read.csv('anuranOccOrig.csv') %>% select(-X)
# count occurrences pulled down ----
occ_count<-anuran_occ %>% 
  group_by(species) %>%
  mutate(n_step5=n()) %>%
  distinct(decimalLatitude, decimalLongitude, year, .keep_all = TRUE) %>%
  mutate(n_step6.1=n()) %>%
  slice(1) %>% select(species, n_step5, n_step6.1)

occ_count %>% arrange(n_step6.1)

#read in each dataset and write out by species ----
main_dir<-'G:/Shared drives/Anuran Biodiversity Project/SDMs/Occurrence Data/Download & Filtering/'
AnuranSp_all<-grep(list.files(main_dir), pattern=".doc", invert=T, value=T)
anuran.taxa.table

# identify the folders that have HerpMapper data in it
hmfolders<-list.files(main_dir, pattern="HerpMapper", 
                      recursive=T, include.dirs = T, full.names=T)
#bind all the HerpMapper data together (should be v3)
library(data.table)
hm_occs<-rbindlist(lapply(hmfolders,
                  function(x){read.csv(list.files(x, 'v3.csv', full.names = T))}))
#clean up the dataset
hm_occs<-hm_occs %>% 
  mutate(source="HM",
         #improve the data column
         DateF=case_when(substr(Date,1,4)=="0000"~as.character(Sys.Date()),
                         substr(Date,9,10)=="00" ~ paste0(substr(Date,1,8),
                                                          as.numeric(substr(Date,9,10))+1),
                         T~Date),
         Date=as.Date(DateF), #ensure Date column is in Date format
         year=year(Date)) %>%  #pull year out from Date column
  select(-DateF) %>%
  rowwise() %>%
  #add species and infraspecificEpithet to match gbif
  mutate(species=paste(scan(text=Taxon, what=" ", n=2, quiet=T), collapse=" "), #paste the species
         infraspecificEpithet=scan(text=Taxon, what=" ", quiet=T)[3]) %>%
         #tax=paste(substr(species,1,1), gsub("[A-z ]*\\ ","", species),sep="."),)%>% #paste the subspecies)
  left_join(anuran.taxa.table) #join for later sorting
gbif_occs<- anuran_occ %>%
  rename("Latitude"="decimalLatitude","Longitude"="decimalLongitude") %>% #match HM dataframe
  mutate(source="GBIF") %>%
  left_join(anuran.taxa.table) #join for later sorting

#count the number of occurrences now in this dataset
fin_counts<-bind_rows(gbif_occs, hm_occs) %>%
  distinct(final.taxa, Longitude, Latitude, year, .keep_all=T) %>%
  count(final.taxa) %>%
  filter(!(final.taxa %in% sp_exclude))

fin_counts %>% arrange(n) #likely should exclude < 10 
write.csv(fin_counts, 'rcs_results/occurrence_counts_used.csv')
#checking that taxa are congruent and which to use for selecting on
bind_rows(gbif_occs, hm_occs) %>% group_by(species, source) %>% slice(1) %>% select(final.taxa,tax,species, family, Longitude, Latitude, year, source,) %>% View()

#for loop to write out individual taxa
for(k in fin_counts$final.taxa){
  tax_ids<-anuran.taxa.table[anuran.taxa.table$final.taxa==k,] #in case there are multiple species
  gbif_k<-gbif_occs %>% filter(species %in% tax_ids$species) #pull all gbif of that species
  hm_k<-hm_occs %>% filter(species %in% tax_ids$species) #pull all hm of that species
  all_k_occ<-bind_rows(gbif_k, hm_k) %>% #bind gbif & hm together
    distinct(Longitude, Latitude, year, .keep_all=T) #remove duplicate location/year data
  print(paste(tax_ids$species)) #print multiple species (as a check)
  write.csv(all_k_occ, 
            paste0('data/occ_data/', k, '_', format(Sys.Date(), "%Y%m%d"), '.csv'),
            row.names = F)
}

unique(hm_occs$species)
unique(gbif_occs$species)
