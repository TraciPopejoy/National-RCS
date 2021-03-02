library(rgbif); library(tidyverse); library(sf)
user="tracidubose" # GBIF user name
pwd="anuran4eva" # GBIF password
email="tracipdubose@gmail.com" # your email
#2427616  #Acris crepitans
# 952 #Anura
gbif_download = occ_download(
  pred("taxonKey",  2427616), 
  pred("hasCoordinate", TRUE), #step 2
  pred("hasGeospatialIssue", FALSE), #partial step 3
  pred_not(pred("issue", "IDENTIFIED_DATE_UNLIKELY")), #step 3
  pred_not(pred("issue", "RECORDED_DATE_MISMATCH")), #step 3
  pred_notnull("country"), #step 4
  pred("continent","North America"), #attempt at step 5
  pred_lte("eventDate", "2019-12-25"),
  format = "SIMPLE_CSV",
  user=user,pwd=pwd,email=email
)
occ_download_wait(gbif_download)
occ_download_get(gbif_download[1], path = "./data")
unzip("./data/0204573-200613084148143.zip",
      exdir="data")
acrep_occ<-read_tsv("./data/0204573-200613084148143.csv")


# old occurrence loss information
# all species info here: https://docs.google.com/spreadsheets/d/1D4QoZGxqvGtTrCSsDcZIsS--QHWQA__eyfiiw21sfE8/edit#gid=1759593941
# Start:30787	
# after step 2: 11728
# after step 3: 11515
# after step 4: 11515
# after step 5: 11489
# after step 6: 2876

#theoretically, we should have same number of occurrences at step 5
nrow(acrep_occ) # more because of september 2020 to dec 2020?
library(lubridate)
acrep_occ %>%
  filter(eventDate <= ymd("2020-09-05")) %>% nrow() # apparently not
#complete step 6
acrep_occ %>%
  distinct(decimalLatitude, decimalLongitude, year, .keep_all = TRUE) %>%
  nrow() #still missing 500 points

#taxonomic movements?

#map where the points occur ----
library(maps)
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
acrep_sf<-acrep_occ %>% st_as_sf(coords=c("decimalLongitude", "decimalLatitude"),
                        crs=4326)
ggplot()+
  geom_sf(data=states)+
  geom_sf(data=acrep_sf, alpha=0.2)+
  coord_sf(xlim=st_bbox(acrep_sf)[c(1,3)],
           ylim=st_bbox(acrep_sf)[c(2,4)])


#ISSUES PRESENT ----
tibble(iss=unique(acrep_occ$issue)) %>% 
  separate(iss, into=paste0("x_", 1:6), sep=";") %>%
  mutate(test="test") %>%
  pivot_longer(-test) %>%
  filter(!is.na(value)) %>%
  group_by(value) %>% slice(1) %>% select(value)

acrep_ndup<-acrepitans %>%
  arrange(desc(year)) %>%
  filter(!duplicated(decimalLatitude, decimalLongitude, year))
nrow(old_acrep);nrow(acrep_ndup)

nrow(acrep)

ranas<-read.csv('data/occ_data_used/Rana muscosa.csv') %>%
  bind_rows(read.csv('data/occ_data_used/Rana pretiosa.csv'))
ranas %>%
  ggplot()+
  geom_histogram(aes(x=year, fill=final.taxa), alpha=0.5) +
  scale_x_continuous(limits=c(1870,2020))+
  facet_grid(~final.taxa, scales='free_y')+
  scale_fill_discrete(guide=F)
