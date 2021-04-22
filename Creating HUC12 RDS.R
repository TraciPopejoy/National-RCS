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
  st_transform(st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84")) %>%
  rename(WSAREA='areasqkm') %>% #rename for convenient fn use later
  filter(!(states %in% c("CN", "MX", "HI", "VI","PR", "AS","GU","CM", "AK", "AK,CN")),
         hutype != 'F')
usa_l48<-st_transform(usa_l48, st_crs(huc12_raw)) #transform to have the same crs
 
# Subset the HydroBASINS used in this analysis (outside contiguous US) ----
hydro12_na_raw<-st_read(dsn=PATH_hydrosheds, layer='hybas_na_lev12_v1c', quiet = T) 
st_crs(hydro12_na_raw)==st_crs(huc12_raw) #check the crs to be the same
ov_ws<-st_intersects(huc12_raw, hydro12_na_raw)

hydro12_na<-hydro12_na_raw[-unique(unlist(ov_ws)),] %>%
  rename(WSAREA='SUB_AREA') #rename for convenient fn use later

# create polygons to bridge area between hydrobasins and hucs ----
# identify polygons that were removed but touch the kept hydro12_na (to reduce shapefile size)
touching<-st_touches(hydro12_na, hydro12_na_raw[unique(unlist(ov_ws)),], sparse=F)
beep<-which(touching, arr.ind=T) #which are touching
#only keep the ones in the raw dataset that were excluded because of hucs
xeep<-hydro12_na_raw[unique(unlist(ov_ws)),] 

ggplot()+#geom_sf(data=hydro12, color='red', fill='red')+
  geom_sf(data=xeep[beep[,2],])

#identify the polygons hucs that were excluded (xeep) and touch the ones that were kept (beep)
overlapped_poly<-xeep[beep[,2],] %>% pull(HYBAS_ID) %>% unique()

rm(ov_ws, xeep, beep, main, US, touching, df) # removing large files I don't need anymore

#keep the polygons that overlap with the edge of the hucs
# (called some_trouble because I need to remove overlapped areas)
some_trouble<-hydro12_na_raw %>%
  filter(HYBAS_ID %in% overlapped_poly) 
# which polygons overlap the hucs?
edge_hydro<-st_intersects(some_trouble, (huc12_raw))

str(edge_hydro)

ggplot()+
  geom_sf(data=some_trouble[426:438,], fill='yellow')+
  geom_sf(data=huc12_raw[unique(unlist(edge_hydro[426:438]))[3],], aes(fill=huc12))

# doing this for each intersection to avoid st_union(huc12), which is a monster to run
library(data.table)
fDT2 <- function(u,n){
  #use rbindlinst to avoid memory slowdown
  return(rbindlist(
    lapply(u:n, #start and stop values for each row of some_trouble
           function(x){
             if(length(edge_hydro[[x]])!=1){ # don't need to do st_union if there is only one polygon
               if(max(huc12_raw[edge_hydro[[x]],]$WSAREA) > 10000){ #large polygons are difficult for st_union
                 print(paste(x, 'had a big polygon')) #identify rows that have a big polygon
                 #remove rows that are a large polygon (generally removed later from huc12)
                 ig.me<-which(huc12_raw[edge_hydro[[x]],]$WSAREA==max(huc12_raw[edge_hydro[[x]],]$WSAREA)) 
                 #identify the area outside the huc12 polygon and shift it to individual polygons
                 st_difference(some_trouble[x,], st_union(huc12_raw[edge_hydro[[x]][-ig.me],])) %>%
                   st_cast("POLYGON")}else{
                     #if no large polygons, just union all huc12 polygons
                     st_difference(some_trouble[x,], st_union(huc12_raw[edge_hydro[[x]],])) %>%
                       st_cast("POLYGON")}}else{
                         #if only one polygon, don't need to st_union at all
                         st_difference(some_trouble[x,], huc12_raw[edge_hydro[[x]],]) %>%
                           st_cast("POLYGON")}}),
    fill=T))
}

# EXAMPLE OF FUNCTION ABOVE 
# when just one polygon exists
p.onep<-st_difference(some_trouble[x,], huc12_raw[edge_hydro[[x]],]) %>%
  ggplot()+
  geom_sf(data=huc12_raw[edge_hydro[[x]],], aes(fill=huc12), alpha=0.3)+
  geom_sf(data=some_trouble[x,], aes(fill='old hydrobasin'), alpha=0.3)+
  geom_sf(aes(fill='kept'))+
  ggtitle('Just one polygon')
# removing a great lake polygon  
p.gl<-st_difference(some_trouble[321,], 
              st_union(huc12_raw[edge_hydro[[321]][-which(huc12_raw[edge_hydro[[x]],]$WSAREA==max(huc12_raw[edge_hydro[[x]],]$WSAREA)],]))) %>%
  ggplot()+
  geom_sf(data=huc12_raw[edge_hydro[[321]],], aes(fill=huc12), alpha=0.3)+
  geom_sf(data=some_trouble[321,], aes(fill='old hydrobasin'), alpha=0.3)+
  geom_sf(aes(fill='kept'))+
  ggtitle('Ignoring great lake polygon')
#when multiple polygons exist
p.multip<-st_difference(some_trouble[1,], st_union(huc12_raw[edge_hydro[[1]],])) %>%
  ggplot()+
  geom_sf(data=huc12_raw[edge_hydro[[1]],], aes(fill=huc12), alpha=0.3)+
  geom_sf(data=some_trouble[1,], aes(fill='old hydrobasin'), alpha=0.3)+
  geom_sf(aes(fill='kept')) +
  ggtitle('Multiple polygons')+
  scale_fill_viridis_d()+
  theme_bw()
plot_grid(p.onep, p.gl, p.multip, nrow=1)

#actually run the function
seq(from=1, to=nrow(some_trouble), nrow(some_trouble)/3)
polys1<-fDT2(1,351)
polys2<-fDT2(352,701)
polys3<-fDT2(702, nrow(some_trouble)) #this one takes a while because of the great lakes

saveRDS(rbind(polys1, polys2, polys3) %>% st_as_sf(), '/home/tracidubose/rcs_newpoly_20210323.rds')

all_new_poly<-rbind(polys1, polys2, polys3) %>% st_as_sf()

#take a look at the results
ggplot()+geom_sf(data=all_new_poly, alpha=0.3)+
  geom_sf(data=huc12_raw[unique(unlist(edge_hydro)),], color='red', alpha=0.3)


# Finish removing large hucs from hucs shapefile (don't want these in here at all)
gl_artifacts<-huc12_raw %>% 
  filter(hutype == "F", #pull out the 'frontal' hucs
         grepl('Lake Superior', name)|grepl('Lake Huron', name)|
           grepl('Lake Michigan', name)|grepl('Lake Erie', name)|
           grepl('Lake Ontario', name), #find hucs with these names
         shape_Length > 20) %>% #keep anything exceptionally large
  pull(huc12)
# I found some more large hucs I want to remove
strange_large_hucs<-huc12_raw %>% 
  filter(grepl('CN', states)) %>% 
  arrange(desc(shape_Length)) %>%
  slice(1:2) %>% pull(huc12) 

huc12<-huc12_raw %>%
  filter(!(name %in% c('Lake Huron', 'Lake Erie', 'Lake Superior', #remove great lakes
                       'Lake Ontario', 'Lake Michigan')),
         !(huc12 %in% c(as.character(gl_artifacts),  #remove strange large hucs
                        as.character(strange_large_hucs))),
         !grepl('AK', states), #remove alaskan hucs
         !is.na(states)) %>% 
  dplyr::select(-noncontributingareaacres,-noncontributingareasqkm) #removing bad column names


# check that new polygons don't overlap huc12 or itself ----
# first check against itself
suspicious<-st_overlaps(all_new_poly, sparse=F)
take_a_look<-which(suspicious, arr.ind=T) # these polygons overlap
#after investigating, they just touch. I'm going to leave them
st_intersection(all_new_poly[unique(c(take_a_look[,1],
                                      take_a_look[,2])),]) %>% 
  filter(n.overlaps!=1) %>% st_area() %>%
  round(3)

#check new polys don't overlap with hydro12
check1<-st_overlaps(all_new_poly, hydro12_na, sparse=F)
(peek1<-which(check1, arr.ind = T))
st_intersection(all_new_poly[peek1[,1],],
                hydro12_na[peek1[,2],]) %>% 
  st_area() %>%
  round(3) #so very small

#check new polys don't overlap with huc12 !very long to run
check2<-st_overlaps(all_new_poly, huc12, sparse=F)
(peek2<-which(check2, arr.ind = T))
huc_areas<-st_intersection(all_new_poly[peek2[,1],],
                huc12[peek2[,2],]) %>% 
  st_area() %>%
  round(3)
ggplot()+
  geom_histogram(aes(x=as.numeric(huc_areas)), bins=20)+
  scale_x_continuous(expression('Overlapped Area m'^2))+
  scale_y_log10()+theme_bw()
ggsave('/home/tracidubose/RCS_Results/sup_fig/Overlapped_Areas_Defaul_poly.jpg',
       width=6, height=4)
units(huc_areas)
#i'm ok with this overlap.


#to do: 
#add area calculation to all_new_poly
names(all_new_poly)
all_new_poly_adj<-all_new_poly %>%
  select(HYBAS_ID:SORT) %>%
    mutate(oSUB_AREA=SUB_AREA,
           SUB_AREA=as.numeric(st_area(.))/1000000) %>%
    group_by(HYBAS_ID) %>%
    mutate(n=n()) %>% 
    #select(-NEXT_DOWN, -NEXT_SINK, -MAIN_BAS, -PFAF_ID) %>% 
    rowid_to_column()  %>%
    ungroup() %>% 
    mutate(HYBAS_ID=case_when(n==1 ~ HYBAS_ID,
                               n != 1 ~ as.numeric(paste(HYBAS_ID, rowid, sep='.')))) %>%
  select(-rowid, oSUB_AREA, -n) %>%
  rename(WSAREA='SUB_AREA')

all_new_poly_adj %>% group_by(HYBAS_ID) %>% count() %>% arrange(desc(n)) #check unique ids
#add arctic hydrobasins
hydro12_ar<-st_read(dsn=PATH_hydrosheds, layer='hybas_ar_lev12_v1c', quiet = T) %>%
  rename(WSAREA='SUB_AREA')
#NEED TO MAKE SURE UNIQUE HYBAS_ID
#rbind hydro12, all_new_poly, and arctic hydrobasing
names(hydro12_na)
names(all_new_poly_adj)
names(hydro12_ar)
hydro12<-rbind(hydro12_ar, hydro12_na, all_new_poly_adj %>% select(-oSUB_AREA))

#some weird hydrobasins snuck in, removing them here

peep<-st_intersects(hydro12, usa_l48 %>% st_buffer(-.3), sparse=F)
hydro12<-hydro12[-which(peep, arr.ind = T)[,1],]

#check these are the correct files
ggplot()+
  #geom_sf(data=huc12, fill='forestgreen', color='forestgreen')+
  geom_sf(data=hydro12, fill='navy', color='navy')+
  geom_sf(data=usa_l48, color='yellow', fill=NA)

#saving copies of the huc12 polygons and USA polygon with the NAD83 projection ----
saveRDS(huc12, '/home/tracidubose/huc12_wgs84_sf_02010323.rds')
saveRDS(hydro12, '/home/tracidubose/hydrobasins_wgs84_sf_20210323.rds')

# Supplemental Fig 1 ----
library(maps)
crs.geo <- st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84")
crs.albers <- st_crs("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs")#2.12kb
dont_plot<-st_as_sf(maps::map('world', c('cuba', 'haiti', 'dominican republic',
                                   'puerto rico','jamaica', 'bahamas', 'dominica',
                                   'anguilla','guadeloupe',
                                   'st lucia', 'cayman islands'), plot = FALSE, fill = TRUE)) %>% 
  st_transform(crs.albers) %>% st_buffer(200) %>% 
  st_transform(crs.geo) %>%
  st_intersection(hydro12) %>%  pull(HYBAS_ID)

spat_ext<-ggplot()+
  geom_sf(data=hydro12[!(hydro12$HYBAS_ID %in% dont_plot),],
          aes(fill='HydroBASINS', color='HydroBASINS'), alpha=0.5)+
  geom_sf(data=huc12,  aes(fill='HUC12', color='HUC12'), alpha=0.5)+
  scale_color_manual(values=c('lightgrey', 'darkgrey'), 
                     aesthetics = c('fill', 'color'),
                     guide=F)+
  theme_bw()+
  theme(axis.text.y=element_text(size=6),
        axis.text.x=element_text(size=6))
examp<-ggplot()+
  geom_sf(data=hydro12, aes(fill='HydroBASINS', color='HydroBASINS'), alpha=0.7)+
  geom_sf(data=all_new_poly, aes(fill='New poly', color='New poly'), alpha=0.3)+
  geom_sf(data=huc12, aes(fill='HUC12', color='HUC12'), alpha=0.7)+
  coord_sf(xlim=c(-97.5, -95.5),
           ylim=c(48, 50))+
  scale_color_manual(values=c('lightgrey','darkgrey','red'), 
                     aesthetics = c('fill', 'color'),
                     guide=F)+
  theme_bw()+
  theme(legend.position='bottom',
        axis.text.y=element_text(size=8),
        axis.text.x=element_text(size=8, angle=20, hjust=.9))
library('cowplot')
plot_grid(spat_ext, examp, labels="AUTO")
ggsave('/home/tracidubose/RCS_Results/sup_fig/supfigspatial20210323.jpeg',
       height=5, width=7) 
ggsave('/home/tracidubose/RCS_Results/sup_fig/Fig2_test_map.jpeg', 
       spat_ext, width=2, height=2)
