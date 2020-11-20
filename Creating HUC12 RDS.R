#create the two huc12 files needed
#need large memory to do this (use HPC available on ARC)
#need to save one as crs nad83/EPSG: 4269 and save one as wgs84
#need to only keep lower 48. 

library(maps); library(tidyverse); library(sf)

# Set external data paths ----
PATH_HUC12 <- '/home/tracidubose/WBD_National_GDB.gdb'
PATH_hydrosheds<-'/home/tracidubose/HydroSHEDS/'

# Create the map of the contiguous US ----
US <- map_data("state") #converts state map to dataframe for mapping with ggplot2
usa <- maps::map("usa") #gets map for subsetting entries. Map of the lower 48 obtained through Package ‘maps’
#version 3.3.0 (Becker et al. 2018) in Program R (R Core Team 2016);
main <- map_data("usa", region=c("main")) #subsets map points to only mainland, no islands

df <- data.frame(main$long, main$lat) #dataframe of longitude and latitude 
usa_l48 <- st_as_sf(df, coords=c("main.long","main.lat"),
                    crs=st_crs("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84")) %>%
  st_combine(.) %>%
  st_cast(., "POLYGON")
# So now, usa_l48 is a map of the contiguous US

# Load, transform, and filter the HUC12 dataset ----
huc12_raw <- st_read(dsn = PATH_HUC12, layer = "WBDHU12") %>%
  st_transform(st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84"))
usa_l48<-st_transform(usa_l48, st_crs(huc12_raw)) #transform to have the same crs
 
gl_artifacts<-huc12_raw %>% 
  filter(hutype == "F",
         grepl('Lake Superior', name)|grepl('Lake Huron', name)|
           grepl('Lake Michigan', name)|grepl('Lake Erie', name)|grepl('Lake Ontario', name),
         shape_Length > 20) %>%
  pull(huc12)
strange_large_hucs<-huc12_raw %>% filter(grepl('CN', states)) %>%
  arrange(desc(shape_Length)) %>%
  slice(1:2) %>% pull(huc12) 

huc12<-huc12_raw %>%
  rename(WSAREA='areasqkm') %>% #rename for convenient fn use later
  filter(!(states %in% c("CN", "MX", "HI", "VI","PR", "AS","GU","CM")),
         !(name %in% c('Lake Huron', 'Lake Erie', 'Lake Superior',
                       'Lake Ontario', 'Lake Michigan')),
         !(huc12 %in% c(as.character(gl_artifacts), 
                        as.character(strange_large_hucs))),
         !grepl('AK', states),
         !is.na(states)) %>%
  dplyr::select(-noncontributingareaacres,-noncontributingareasqkm) #removing bad column names


# Subset the HydroBASINS used in this analysis (outside contiguous US) ----
hydro12_raw<-st_read(dsn=PATH_hydrosheds, layer='hybas_na_lev12_v1c', quiet = T) 
st_crs(hydro12_raw)==st_crs(huc12)
ov_ws<-st_intersects(huc12, hydro12_raw)

hydro12<-hydro12_raw[-unique(unlist(ov_ws)),] %>%
  rename(WSAREA='SUB_AREA') #rename for convenient fn use later

#saving copies of the huc12 polygons and USA polygon with the NAD83 projection ----
saveRDS(usa_l48, '/home/tracidubose/usal48_nad83_sf.rds')
saveRDS(huc12, '/home/tracidubose/huc12_wgs84_sf.rds')
saveRDS(hydro12, '/home/tracidubose/hydrobasins_wgs84_sf.rds')
