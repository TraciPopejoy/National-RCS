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
usa_l48<-st_transform(usa_l48, st_crs(huc12))

huc12_alb<-st_transform(huc12, 
                       st_crs("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs"))
usa_l48_alb<-st_transform(usa_l48, st_crs(huc12_alb))

huc12_filt_cont <- huc12_alb %>%
  filter(st_contains(usa_l48_alb, ., sparse = FALSE))

pdf('/home/tracidubose/huc12_filt_cont.pdf')
plot(st_geometry(huc12_filt_cont))
dev.off()

huc12_crop<-st_crop(huc12_alb, st_bbox(usa_l48_alb))
huc12_crop_old<-huc12_crop %>% filter(states!="MX", states!="CN",
                                      !(name %in% c("Lake Superior","Lake Erie","Lake Ontario","Lake Huron","Lake Michigan")))


# not working, can't get sf to update
#library(devtools)
#install_github("r-spatial/sf") #need sf version > 8.0 (works at 9.5) for st_filter
#huc12_snap<-sf::st_filter(huc12_alb,
#                     st_as_sf(usa_l48_albs))


#saving copies of the huc12 polygons and USA polygon with the NAD83 projection
saveRDS(usa_l48, '/home/tracidubose/usal48_nad83_sf.rds')
saveRDS(huc12_crop, '/home/tracidubose/huc12_nad83_sf_crop.rds')
saveRDS(huc12_filt_cont, '/home/tracidubose/huc12_nad83_sf_mask.rds')
huc12_geo<-st_transform(huc12, st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84"))
saveRDS(huc12_geo, '/home/tracidubose/huc12_wgs84_sf.rds')
saveRDS(huc12_alb, '/home/tracidubose/huc12_albers_sf.rds')

#look to see if mask actually worked well

hc12mex_mask<-huc12_filt_cont %>% filter(grepl('MX', states))
all_hc12mex<-huc12_alb %>% filter(grepl('MX', states))
hc12mex_crop<-huc12_crop_old %>% filter(grepl('MX', states))

st_bbox(hc12mex)

ggplot()+
  geom_sf(data=all_hc12mex, fill="blue", alpha=0.3)+
  geom_sf(data=hc12mex_crop, fill="yellow", alpha=0.5) +
  geom_sf(data=hc12mex_mask, fill="green", alpha=0.5) +
  geom_sf(data=usa_l48_alb, color="white", size=2, fill=NA)+
  coord_sf(xlim=st_bbox(hc12mex[1:2,])[c(1,3)],
           ylim=st_bbox(hc12mex[1:2,])[c(2,4)])
  