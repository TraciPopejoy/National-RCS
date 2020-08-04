#create the two huc12 files needed
#need large memory to do this (use HPC available on ARC)
#need to save one as crs nad83/EPSG: 4269 and save one as wgs84
#need to only keep lower 48. 

library(maps); library(tidyverse); library(sf)

# Create the map of the contiguous US
US <- map_data("state")#converts state map to dataframe for mapping with ggplot2
usa <- maps::map("usa") #gets map for subsetting entries. Map of the lower 48 obtained through Package ‘maps’
#version 3.3.0 (Becker et al. 2018) in Program R (R Core Team 2016);
main <- map_data("usa", region=c("main"))#subsets map points to only mainland, no islands

df <- data.frame(main$long, main$lat) #dataframe of longitude and latitude 
usa_l48 <- st_as_sf(df, coords=c("main.long","main.lat"),
              crs=st_crs("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84")) %>%
  st_combine(.) %>%
  st_cast(., "POLYGON")

huc12 <- st_read(dsn = '/home/tracidubose/WBD_National_GDB.gdb', layer = "WBDHU12")
usa_l48<-st_transform(p, st_crs(huc12))
huc12_crop<-st_crop(huc12, st_bbox(p))


huc12_crop<-huc12_crop %>% filter(states!="MX", states!="CN",
                      !(name %in% c("Lake Superior","Lake Erie","Lake Ontario","Lake Huron","Lake Michigan")))
huc12_crop %>% arrange(desc(areasqkm)) %>% slice(1:5)

#pdf('/home/tracidubose/huc12_included.pdf')
#plot(st_geometry(huc12_crop))
#dev.off()

#saving copies of the huc12 polygons and USA polygon with the NAD83 projection
saveRDS(usa_l48, '/home/tracidubose/usal48_nad83_sf.rds')
saveRDS(huc12_crop, '/home/tracidubose/huc12_nad83_sf.rds')
huc12_geo<-st_transform(huc, st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84"))
saveRDS(huc12_geo, '/home/tracidubose/huc12_wgs84_sf.rds')

#### investigating spatial extent, sf, and mask ####
huc2 <- st_read(dsn = 'C:/Users/Owner/Documents/GISfile/WBD_National_GDB/WBD_National_GDB.gdb', 
                 layer = "WBDHU2")
huc2_alb<-st_transform(huc2, 
                       st_crs("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs"))
p<-st_transform(p, st_crs(huc2_alb))

huc2_filt_cont <- huc2_alb %>% 
  st_buffer(dist=0) %>%
  filter(st_contains(uusa_l48_albs, ., sparse = FALSE))

huc2_crop<-st_crop(st_buffer(huc2_alb, dist=0),
                   st_bbox(p))
huc2_snap<-st_filter(st_buffer(huc2_alb, dist=0),
                     st_as_sf(usa_l48_albs))

ggplot()+
  geom_sf(data=p)+
  geom_sf(data=huc2_crop, alpha=0.5, fill="goldenrod")+
  geom_sf(data=huc2_snap, alpha=0.3, fill="orange")