#create the two huc12 files needed
#need large memory to do this (use HPC available on ARC)
#need to save one as crs nad83/EPSG: 4269 and save one as wgs84
#need to only keep lower 48. 

library(maps); library(tidyverse); library(sf)

US <- map_data("state")#converts state map to dataframe for mapping with ggplot2
usa <- maps::map("usa") #gets map for subsetting entries. Map of the lower 48 obtained through Package ‘maps’
#version 3.3.0 (Becker et al. 2018) in Program R (R Core Team 2016);
main <- map_data("usa", region=c("main"))#subsets map points to only mainland, no islands

#converts from dataframe to spatial polygon
df <- data.frame(main$long, main$lat) #dataframe of longitude and latitude 
p <- st_as_sf(df, coords=c("main.long","main.lat"),
              crs=st_crs("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84")) %>%
  st_combine(.) %>%
  st_cast(., "POLYGON")

huc12 <- st_read(dsn = '/home/tracidubose/WBD_National_GDB.gdb', layer = "WBDHU12")
p<-st_transform(p, st_crs(huc12))
huc12_crop<-st_crop(huc12, st_bbox(p))


huc12_crop<-huc12_crop %>% filter(states!="MX", states!="CN",
                      !(name %in% c("Lake Superior","Lake Erie","Lake Ontario","Lake Huron","Lake Michigan")))
huc12_crop %>% arrange(desc(areasqkm)) %>% slice(1:5)

#pdf('/home/tracidubose/huc12_included.pdf')
#plot(st_geometry(huc12_crop))
#dev.off()

#saving a copy with the NAD83 projection
saveRDS(huc12_crop, '/home/tracidubose/huc12_nad83_sf.rds')
huc12_geo<-st_transform(huc, st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84"))
saveRDS(huc12_geo, '/home/tracidubose/huc12_wgs84_sf.rds')
