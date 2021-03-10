library(rgbif); library(tidyverse); library(sf)
library(CoordinateCleaner) #to remove sea coordinates
user="tracidubose" # GBIF user name
pwd="anuran4eva" # GBIF password
email="tracipdubose@gmail.com" # your email

# don't pull down occurrences for species I'm not going to use ----
# excluding those with < 20% range in US
sp_range_wiL48 <-read.csv('rcs_results/sp_range_within_US.csv') %>% 
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

# old occurrence loss information ----
# all species info here: https://docs.google.com/spreadsheets/d/1D4QoZGxqvGtTrCSsDcZIsS--QHWQA__eyfiiw21sfE8/edit#gid=1759593941
orig_filt<-read.csv("C:/Users/Owner/Downloads/Occurrence Record Tracksheet - ALL SPECIES.csv")
names(orig_filt) <- orig_filt[1,]
orig_filt<-orig_filt[-c(1,110:114),-c(9:12)]

orig_filt<-orig_filt %>% 
  mutate(tax=paste0(substr(Species,1,1),'.',gsub(".*? ","",Species))) %>%
  left_join(anuran.taxa.table, by='tax') %>% 
  filter(!is.na(final.taxa),
         !(species %in% sp_exclude)) %>%
  select(Species, final.taxa, `5`,`6`) %>%
  arrange(desc(`5`)) %>%
  filter(!duplicated(final.taxa))

for(i in 1:nrow(orig_filt)){orig_filt$roll_occ[i]<-sum(orig_filt[1:i,3])}
head(orig_filt)

# get gbif species keys ----
anuran.gbif.taxa<-NULL
for(u in 1:nrow(sp_range_wiL48)){
  anuran.gbif.taxa1<-name_backbone(name=sp_range_wiL48[u,1], 
                                   rank='species', kingdom='animals')
  anuran.gbif.taxa<-bind_rows(anuran.gbif.taxa1, anuran.gbif.taxa)
}
anuran.gbif.taxa<-anuran.gbif.taxa %>%
  filter(!(species %in% sp_exclude), #remove species not interested in 
         !duplicated(speciesKey)) %>% #remove duplicated species to not collect twice
  left_join(orig_filt, by=c('canonicalName'='final.taxa')) %>% #order to optimize download
  arrange(desc(`5`)) %>%
  filter(!duplicated(speciesKey))

View(anuran.gbif.taxa %>% select(species, `5`,roll_occ) %>% left_join(occ_count))

# set up a four loop to pull down occ (> 100k occ) ----
# automate pulling from species list
ind_mat<-matrix(c(1,3,4,11,12,24,25,92), ncol=2, byrow = T)
anuran_occ<-NULL #where all occ kept eventually
gbif_codes<- NULL #keep the download codes

for(j in 1:5){
  gbif_download = occ_download(
    #step 1 - get occ for a species
    pred_in("taxonKey",  anuran.gbif.taxa$speciesKey[ind_mat[j,1]:ind_mat[j,2]]), 
    pred("hasCoordinate", TRUE), #step 2
    pred("hasGeospatialIssue", FALSE), #partial step 3
    pred_not(pred("issue", "IDENTIFIED_DATE_UNLIKELY")), #step 3
    pred_not(pred("issue", "RECORDED_DATE_MISMATCH")), #step 3
    pred_notnull("country"), #step 4
    pred_in("country",c("US","MX","CA")), #attempt at step 5
    pred_lte("eventDate", "2019-12-25"),
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

doublecheck<-NULL
for(k in gbif_codes){
  assign(paste0('occ_',grep(k, gbif_codes)),
         read_tsv(paste0("./data/", k,".csv"),
                 col_types="cccccccccccccccccccicnnnnnnnnTnnnccccccccTcccccTcc"))
  doublecheck<-bind_rows(doublecheck, get(paste0('occ_',grep(k, gbif_codes))))
}

length(c(unique(occ_1$species), #three lithobates
unique(occ_2$species), #four
unique(occ_3$species),
unique(occ_4$species),
unique(occ_5$species)))

View(anuran_occ)

# remove occ with NA
na_occ<-doublecheck %>% filter(is.na(decimalLatitude))

# remove occurrences within the ocean 
anuran_occ_orig<-anuran_occ
doublecheck <- doublecheck %>%
  filter(!is.na(decimalLatitude)) %>%
  cc_sea(lon = "decimalLongitude",  lat = "decimalLatitude")

# count occurrences pulled down ----
occ_count<-doublecheck %>% 
  group_by(species) %>%
  mutate(n_step5=n()) %>%
  distinct(decimalLatitude, decimalLongitude, year, .keep_all = TRUE) %>%
  mutate(n_step6.1=n()) %>%
  slice(1) %>% select(species, n_step5, n_step6.1)



anuran.taxa.table<-read.csv('./data/AnuranTaxRef_20200708.csv') %>% 
  right_join(occ_count, by=c("final.taxa"="species"))

occ_difference<-occ_count %>% left_join(orig_filt, by=c('species'='final.taxa')) %>%
  select(species, `5`, n_step5, `6`, n_step6.1) %>%
  mutate(step5_dif=`5`-n_step5,
         step6_dif=`6`-n_step6.1)

occ_difference %>%
  ggplot()+geom_histogram(aes(x=step6_dif))

occ_difference %>% arrange(step6_dif)
