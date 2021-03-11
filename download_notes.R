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

head(orig_filt)

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
  
for(k in 1:nrow(anuran.gbif.taxa)){
  anuran.gbif.taxa$n_occ[k]<-occ_count(taxonKey=anuran.gbif.taxa$speciesKey[k], 
                                        georeferenced = T)
  anuran.gbif.taxa <-arrange(anuran.gbif.taxa, desc(n_occ))
  anuran.gbif.taxa$roll_occ[k]<-sum(anuran.gbif.taxa$n_occ[1:k])
}

View(anuran.gbif.taxa %>% select(canonicalName, speciesKey, n_occ, roll_occ))

# set up a four loop to pull down occ (> 100k occ) ----
# automate pulling from species list
ind_mat<-matrix(c(1,3,4,11,12,24,25,92), ncol=2, byrow = T)
anuran_occ<-NULL #where all occ kept eventually
gbif_codes<- NULL #keep the download codes

for(j in 1:4){
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
write.table(gbif_codes,
            "GBIF_codes_20210311.txt")
View(anuran_occ)

# remove occ with NA
na_occ<-anuran_occ %>% filter(is.na(decimalLatitude))

# remove occurrences within the ocean 
anuran_occ_orig<-anuran_occ
anuran_occ <- anuran_occ %>%
  filter(!is.na(decimalLatitude)) %>%
  cc_sea(lon = "decimalLongitude",  lat = "decimalLatitude")

# count occurrences pulled down ----
occ_count<-anuran_occ %>% 
  group_by(species) %>%
  mutate(n_step5=n()) %>%
  distinct(decimalLatitude, decimalLongitude, year, .keep_all = TRUE) %>%
  mutate(n_step6.1=n()) %>%
  slice(1) %>% select(species, n_step5, n_step6.1)

occ_difference<-occ_count %>% 
  left_join(orig_filt %>% left_join(anuran.taxa.table, by=c("Species"="species")),
            by=c('species'='final.taxa')) %>%
  select(species, `5`, n_step5, `6`, n_step6.1) %>%
  mutate(step5_dif=`5`-n_step5,
         step6_dif=`6`-n_step6.1)
write_csv(occ_difference,"occ_download_occurrence_differences.csv")
# visualize the difference ----
library(scales)
occ_difference %>%
  ggplot()+
  geom_histogram(aes(x=step6_dif))+
  geom_text(data=data.frame(),
            aes(x=-100, y=10, label="more occ\nnow"))+
  geom_text(data=data.frame(),
            aes(x=10, y=10, label="more occ\npast"))+
  scale_x_continuous("Kept Occurrences", 
                     breaks=c(-3000,-1000,-100,-10,0,10,100,1000,3000),
                     trans=scales::pseudo_log_trans(sigma=.2))+
  theme_bw()+
  theme(axis.text.x = element_text(angle=20, hjust=0.9))

occ_difference %>% ungroup() %>% arrange(step6_dif) %>%
  filter(!is.na(step6_dif)) %>% slice(1:2, (n()-1), n())

#identify species to look at ----
occ_difference %>% arrange(step6_dif)
random_number<-rpois(1,92/5)
rcs_sample<-occ_difference %>% arrange(step6_dif) %>%
  ungroup() %>%
  slice(random_number*1:5)
rcs_sample

#investigating difference in datasets
occ_difference %>% filter(step6_dif < 0) %>% 
  ungroup() %>% sample_n(1)
crs.geo <- st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84")
old_anaxhous<-read.csv('data/occ_data_used/Anaxyrus houstonensis.csv')
new_anaxhous<-anuran_occ %>%
  filter(species == "Anaxyrus houstonensis") %>%
  distinct(decimalLatitude, decimalLongitude, year, .keep_all = TRUE) 

ggplot()+
  geom_sf(data=old_anaxhous %>%
            st_as_sf(coords=c("Longitude","Latitude"), crs = crs.geo),
          color="red", alpha=0.2)+
  geom_sf(data=new_anaxhous%>%
            st_as_sf(coords=c("decimalLongitude","decimalLatitude"), crs = crs.geo),
          color='goldenrod3', alpha=0.2)+
  theme_bw()

old_anaxhous %>% select(Longitude, Latitude, year) %>%
  semi_join(new_anaxhous %>% 
              select(decimalLongitude, decimalLatitude, year) %>%
              rename(Longitude=decimalLongitude, Latitude=decimalLatitude))
# :( none match at all.

crs.albers <- st_crs("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs")#2.12kb
o_ah_buf<-old_anaxhous %>% 
  st_as_sf(coords=c("Longitude","Latitude"), crs = crs.geo) %>%
  st_transform(crs.albers) %>% st_buffer(1) %>%
  st_union()
n_ah_buf<-new_anaxhous %>%
  st_as_sf(coords=c("decimalLongitude","decimalLatitude"), crs = crs.geo) %>%
  st_transform(crs.albers) %>% st_buffer(1) %>%
  st_union()

library(maps)
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE)) %>% 
  st_transform(crs.albers)
ggplot()+
  #geom_sf(data=states[states$ID=='texas',])+
  geom_sf(data=o_ah_buf, fill='red', color='red',alpha=0.3)+
  geom_sf(data=n_ah_buf, fill='goldenrod3', color='goldenrod3', alpha=0.3)

st_area(n_ah_buf)-st_area(o_ah_buf) #change in area

ah_dif_buff<-st_difference(n_ah_buf, o_ah_buf)
ah_dif_buff
ggplot()+geom_sf(data=ah_dif_buff, fill='purple', color='purple')
