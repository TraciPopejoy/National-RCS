# Area of Occupancy Code for the National Fishes Vulnerability Assessment - "Area of Occupancy Code_SCS.R"
# Revised by Sam Silknetter, 03April2020
# edited on 2020Sept5 by TPD

library(sf); library(tidyverse); library(scales) # load necessary libraries.
library(cowplot) #for plotting maps later

# Set data paths -----
#shapefiles used
PATH_iucn_ranges <-'/home/tracidubose/RCS_Anuran_input/ANURA/' #IUCN range maps downloaded 10/2020
PATH_HUC12 <- '/home/tracidubose/huc12_wgs84_sf_02010323.rds'
PATH_hydrosheds<-'/home/tracidubose/hydrobasins_wgs84_sf_20210323.rds'
PATH_l48_polygon <- '/home/tracidubose/usal48_nad83_sf.rds'
PATH_FocalSpecies_OccurrenceData <- "/home/tracidubose/RCS_Anuran_input/Anuran Occurrence Data/"
#output paths for taxa shapefiles
PATH_SHP_HUC12 <- "/home/tracidubose/RCS_Results/20210330_sp_shp_files/huc12_sp/"
PATH_SHP_L48_1km <- "/home/tracidubose/RCS_Results/20210330_sp_shp_files/L48_1km/"
PATH_SHP_HyrdoBasins <- "/home/tracidubose/RCS_Results/20210330_sp_shp_files/NA_hydrobasins/"
PATH_SHP_NA_1km <- "/home/tracidubose/RCS_Results/20210330_sp_shp_files/NA_1km/"

# Set spatial coordinate reference system (CRS).
crs.geo <- st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84")
crs.albers <- st_crs("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs")

# Read in fixed data -----
#IUCN maps for native ranges
iucn_anuran<-read_sf(PATH_iucn_ranges, 'ANURA')  %>% 
  st_transform(crs.albers)
#USGS HUC12 subbasins
huc12<-readRDS(PATH_HUC12) 
st_crs(huc12)==crs.geo
#Lower 48 to split points for prism
usa_l48_alb<-readRDS(PATH_l48_polygon) %>%
  st_transform(., crs=crs.albers)
st_crs(usa_l48_alb)==crs.albers
usa_l48_geo<-st_transform(usa_l48_alb, crs.geo)
#HydroBasins and made polygons to fill gaps
hydro12<-readRDS(PATH_hydrosheds) 
st_crs(hydro12)==crs.geo

#check these are the correct files
#ggplot()+
#  geom_sf(data=huc12, fill='forestgreen', color='forestgreen')+
#  geom_sf(data=hydro12, fill='navy', color='navy')+
#  geom_sf(data=usa_l48_geo, color='yellow', fill=NA)


# Create a list of occurrence data files for all focal species ----- 
tax_reference<-read.csv("/home/tracidubose/RCS_Anuran_input/gbif_anuran_taxa_info.csv")
nrow(tax_reference)
anuran.taxa<-tax_reference$species

find_area<-function(taxa, # which species to find the aoo for
                    spatial.file1, # which sp object to use for filtering/grain sizes
                    anti.spatial.file=F, #do you want everything outside of that spatial.file1?
                    output.path, #where do you want to save the output spatial files?
                    native_range=T, #do you want to exclude points >100km outside the native range
                    watershed=F, #are you using watersheds as grain sizes?
                    rad=1, #if buffering points, what should the buffer radius be?
                    watershed.name.column){
  if(!is.na(spatial.file1)){spatial.file<-get(spatial.file1)}
  #read in sp occurrence data
  geodata<-read.csv(list.files(PATH_FocalSpecies_OccurrenceData, pattern=taxa, full.names = T)) %>%
    dplyr::select(family, genus, species, Longitude, Latitude, year, source, final.taxa)
  
  # Convert from 'geodata' CSV to a sf object. 
  dat_sp <- st_as_sf(geodata, coords=c("Longitude","Latitude"), crs = crs.geo) #Save spatial dataframe with correct projection.
  scientific_name <- as.character(dat_sp$final.taxa[1])
  
  if(native_range==T){
    nat_range <- iucn_anuran %>% filter(binomial %in% taxa) %>%
      st_union() %>% st_buffer(100) 
    dat_sp <- dat_sp %>% st_transform(crs.albers) %>% st_filter(nat_range)
  }
  
  # Calculate area occupied per species. 
  if(watershed==T){
    dat_sp <- st_transform(dat_sp, st_crs(spatial.file)) # ensuring transformations fit between spatial file and points
    area.type<- paste0(spatial.file1, '_WS')
    dat_sp <- st_join(dat_sp, spatial.file)  %>% # Store the watershed name as an attribute of the data.
      rename(watershed.name=all_of(watershed.name.column)) # renaming watershed name column
    N_watersheds <- unique(dat_sp$watershed.name) #identify which watersheds are occupied
    na_pts<-length(is.na(dat_sp$watersed.name))

    finsh_sp <- spatial.file %>%
      rename(watershed.name=all_of(watershed.name.column)) %>%
      filter(watershed.name %in% N_watersheds) # Extract unique HUC12 data for the species.
    area <- sum(finsh_sp$WSAREA, na.rm=T)  # Sum the total area (square kilometer) for each unique HUC12.
    st_write(finsh_sp, dsn = paste0(output.path, geodata$final.taxa[1], '.shp'),
             append=F)
    print(paste0('Number of NA areas =', sum(is.na(finsh_sp$WSAREA))))
    print(paste(length(!is.na(N_watersheds)), class(st_geometry(finsh_sp))))
  }else{
    # Generate 1km point buffers. You must use a projected CRS. For CRS, we use: USA Contiguous albers equal area. 
    dat_sp <- st_transform(dat_sp, crs.albers) #overwrites to the projection 
    
    # Filter to the spatial extent / spatial.file provided
    if(anti.spatial.file==T){
      area.type<- paste0('NOT', spatial.file1,'_',rad,'km')
      dat_sp<-dat_sp %>% filter(!st_contains(spatial.file, ., sparse = FALSE))
    }else{
      area.type<- paste0(spatial.file1,'_',rad,'km')
      dat_sp<-dat_sp %>% filter(st_contains(spatial.file, ., sparse = FALSE))
    }
    # Create 1 km buffer and calculate the total area occupied per species.
    finsh_sp <- st_union(st_buffer(dat_sp, dist = rad))
    if(nrow(dat_sp)==0){
      area<-0
      na_pts<-0
    }else{
      area <- as.numeric(st_area(finsh_sp))
      na_pts<-0
      # Save the spatial dataframes as ESRI Shapefiles. 
      st_write(finsh_sp, dsn = paste0(output.path, geodata$final.taxa[1], '.shp'),
               append=F)}
  }
  AOOs <- data.frame(scientific_name = scientific_name, 
                     area.type = area.type,
                     area.sqkm = area,
                     na.pts=na_pts,
                     stringsAsFactors = F)
  return(AOOs)
}


# for-loop to create 4 different shapefiles for each anuran taxa -----
AOOs<-NULL
#short_run<-c("Lithobates yavapaiensis", "Lithobates onca", "Dryophytes wrightorum")
for(y in anuran.taxa){
  l48_buff_aoo<-find_area(taxa=y, spatial.file1='usa_l48_alb',
                     output.path=PATH_SHP_L48_1km,
                     watershed.name.column = NA)
  huc12_aoo<-find_area(taxa=y, spatial.file1='huc12',
                     output.path=PATH_SHP_HUC12, watershed=T,
                     watershed.name.column = 'huc12')
  hydroNA_aoo<-find_area(taxa=y, spatial.file1='hydro12',
                       output.path=PATH_SHP_HyrdoBasins, watershed=T,
                       watershed.name.column = 'HYBAS_ID')
  NA_buff_aoo<-find_area(taxa=y ,spatial.file1='usa_l48_alb',
                         anti.spatial.file = T,
                       output.path=PATH_SHP_NA_1km, watershed=F,
                       watershed.name.column = NA)
  AOOs<-bind_rows(AOOs, 
                  l48_buff_aoo, NA_buff_aoo,
                  hydroNA_aoo, huc12_aoo)
}
#warnings about long character strings in shapefiles 
#Anaxyrus blanchardi doesn't have any geometries :(

AOOs %>% group_by(scientific_name) %>% summarize(total.area=sum(area.sqkm)) %>%
  filter(total.area < 10)
AOOs_orig<-AOOs
AOOs<-AOOs %>% filter(!(scientific_name %in% c("Acris blanchardi","Lithobates kauffeldi"))) 

for(y in c('Acris blanchardi', 'Lithobates kauffeldi')){
  l48_buff_aoo<-find_area(taxa=y, spatial.file1='usa_l48_alb', native_range = F,
                          output.path=PATH_SHP_L48_1km,
                          watershed.name.column = NA)
  huc12_aoo<-find_area(taxa=y, spatial.file1='huc12', native_range = F,
                       output.path=PATH_SHP_HUC12, watershed=T,
                       watershed.name.column = 'huc12')
  hydroNA_aoo<-find_area(taxa=y, spatial.file1='hydro12', native_range = F,
                         output.path=PATH_SHP_HyrdoBasins, watershed=T,
                         watershed.name.column = 'HYBAS_ID')
  NA_buff_aoo<-find_area(taxa=y ,spatial.file1='usa_l48_alb', native_range = F,
                         anti.spatial.file = T,
                         output.path=PATH_SHP_NA_1km, watershed=F,
                         watershed.name.column = NA)
  AOOs<-bind_rows(AOOs, 
                  l48_buff_aoo, NA_buff_aoo,
                  hydroNA_aoo, huc12_aoo)
}

AOOs %>% group_by(scientific_name) %>% summarize(total.area=sum(area.sqkm)) %>%
  filter(total.area < 10)
AOOs %>% filter(scientific_name %in% c('Acris blanchardi', 'Lithobates kauffeldi'))

head(AOOs)
write.csv(AOOs, paste0('/home/tracidubose/RCS_Results/AOOs_raw_',
                       format(Sys.Date(), "%Y%m%d"),'.csv'))

# Summarize the Area of Occurrence -----
AOOs_summarized<- AOOs %>% 
  #identify which grain type each category is
  mutate(grain_type=case_when(grepl('km', area.type)~'buffer',
                              grepl('WS', area.type)~'watershed'),
         alt_area.sqkm=case_when(area.type=='hydro12_WS'~area.sqkm*0.736, #area reduction for larger hydrobasin 
                                 T~area.sqkm)) %>%
  group_by(scientific_name, grain_type) %>%
  #for each grain type, sum the area to get total area occupied
  dplyr::summarize(total.area.occ.sqkm=sum(alt_area.sqkm)) %>%
  pivot_wider(names_from='grain_type', values_from = 'total.area.occ.sqkm') %>%
  ungroup()%>%
  #add ranks for quick comparison
  mutate(rank_WS=rank(watershed),
         rank_buff=rank(buffer),
         rescale(AOOs$mean_Rank, to = c(0,1))) %>%
  rowwise() %>%
  mutate(mean_rank=mean(c(rank_WS, rank_buff)),
         sd_rank=sd(c(rank_WS, rank_buff)),
         log10_buff_sqkm=log10(buffer),
         log10_WS_sqkm=log10(watershed))

# Write CSV out
write.csv(AOOs_summarized, file = paste0("RCS_Results/AOO HUC12 Output_",
                              format(Sys.Date(), '%Y%m%d'),'.csv'),
          row.names = F)

# CHECK THIS WORKED ------
AOOs %>%  mutate(grain_type=case_when(grepl('km', area.type)~'buffer',
                                      grepl('WS', area.type)~'watershed'))%>%
  filter(grain_type=='buffer') %>%
  pivot_wider(names_from='area.type', values_from = 'area.sqkm') %>%
  mutate(inUSratio=usa_l48_alb_1km/NOTusa_l48_alb_1km) %>%
  filter(inUSratio >1, inUSratio<3) %>%
  arrange(desc(inUSratio))

onespp_Huc12<-read_sf(list.files(PATH_SHP_HUC12, pattern='sicolor.shp', full.names = T))
onespp_L48buf<-read_sf(list.files(PATH_SHP_L48_1km, pattern='sicolor.shp', full.names=T))
onespp_hydro<-read_sf(list.files(PATH_SHP_HyrdoBasins, pattern='sicolor.shp', full.names = T))
onespp_NAbuf<-read_sf(list.files(PATH_SHP_NA_1km, pattern='sicolor.shp', full.names = T))
ggplot()+
  geom_sf(data=usa_l48_alb)+
  geom_sf(data=onespp_hydro, aes(color='HydroBasins',fill='HydroBasins'), alpha=0.2)+
  geom_sf(data=onespp_Huc12, aes(color='HUC 12', fill='HUC 12'),alpha=0.2)+
  geom_sf(data=onespp_NAbuf, aes(color='pts NA',fill='pts NA'), alpha=0.2)+
  geom_sf(data=onespp_L48buf, aes(color='pts L48', fill='pts L48'), alpha=0.2)+
  scale_fill_viridis_d(name="Extent", aesthetics = c('fill', 'color'), end=.8)+
  coord_sf(xlim=st_bbox(onespp_NAbuf)[c(1,3)], 
          ylim=st_bbox(onespp_NAbuf)[c(2,4)])
library(scales)
viridis_pal()(4)
library(tmap)
tmap_mode('view')
tm_shape(onespp_hydro)+tm_polygons(col="#440154FF")+
  tm_shape(onespp_Huc12)+tm_polygons(col="#31688EFF")+
  tm_shape(onespp_NAbuf)+tm_polygons(col="#35B779FF")+
  tm_shape(onespp_L48buf)+tm_polygons(col="#FDE725FF")
  

#supplementary figures -----

geodata<-read.csv(list.files(PATH_FocalSpecies_OccurrenceData, pattern=anuran.taxa[9], full.names = T)) %>%
  dplyr::select(family, genus, species, Longitude, Latitude, year, source, final.taxa)

# Convert from 'geodata' CSV to a sf object. 
dat_sp <- st_as_sf(geodata, coords=c("Longitude","Latitude"), crs = crs.geo)
iucn_anuran %>% filter(binomial %in% anuran.taxa[9]) %>%
  st_union() %>% st_buffer(100) %>%
  ggplot()+
  geom_sf(data=usa_l48_alb)+
  geom_sf(fill='blue', alpha=0.5)+
  geom_sf(data=st_transform(st_union(iucn_anuran[iucn_anuran$binomial==anuran.taxa[9],]), crs.albers),
          fill='grey', alpha=0.8)+
  geom_sf(data=dat_sp, alpha=0.15)+
  ggtitle(anuran.taxa[9])+
  theme_bw()+
  theme(title = element_text(face='italic'))
ggsave('/home/tracidubose/rcs_results/sup_fig/natrange_buf_ex.jpg',
       width=4, height=4)


#compare base sizes of huc12 and hydro12 datasets
huc12_areas_all<-st_area(huc12)
hydro12_areas_all<-st_area(hydro12)
ggplot()+geom_density(aes(x=as.numeric(huc12_areas_all)/1e6, fill="HUC 12"), 
                      alpha=0.5)+
  geom_density(aes(x=as.numeric(hydro12_areas_all)/1e6, fill="HydroBasins"),
               alpha=0.5)+
  scale_x_continuous(expression('Watershed area km'^2),
                trans='log10', labels=scales::breaks_log())+
  scale_colour_viridis_d('Dataset\nOrigin', 
                         aesthetics = c('fill', 'color'),
                         end=.9)+
  theme_bw()
ggsave('/home/tracidubose/rcs_results/sup_fig/WSareasdis.jpg',
       width=4.5, height=2)
t.test(as.numeric(huc12_areas_all), as.numeric(hydro12_areas_all))
mean(as.numeric(huc12_areas_all)) / mean(as.numeric(hydro12_areas_all))

test<-read_sf(PATH_SHP_HUC12, "Ascaphus montanus")

test2<-read_sf(PATH_SHP_HyrdoBasins, "Ascaphus montanus")
ggplot()+geom_sf(data=test, aes(fill='HUC12', color='HUC12'))+
  geom_sf(data=test2, aes(fill='HydroBasins', color="HydroBasins"))+
  scale_colour_viridis_d('Dataset\nOrigin', 
                                   aesthetics = c('fill', 'color'),
                                   end=.9)+
  geom_sf(data=iucn_anuran[iucn_anuran$binomial=="Ascaphus montanus",], alpha=0.5)+
  ggtitle("Dryophytes wrightorum")+
  theme_bw()+
  theme(plot.title = element_text(face='italic'),
        axis.text.y=element_text(size=8),
        axis.text.x=element_text(size=8, angle=20, hjust=.9))
ggsave('/home/tracidubose/RCS_Results/sup_fig/huc example plot iucn range.jpg', width=4.5, height=4)


# Count obs loss ----
counts_final<-NULL
for(k in anuran.taxa){
  geodata<-read.csv(list.files(PATH_FocalSpecies_OccurrenceData, pattern=k, full.names = T)) %>%
    dplyr::select(family, genus, species, Longitude, Latitude, year, source, final.taxa)
  
  # Convert from 'geodata' CSV to a sf object. 
  dat_sp <- st_as_sf(geodata, coords=c("Longitude","Latitude"), crs = crs.geo) #Save spatial dataframe with correct projection.
  scientific_name <- as.character(dat_sp$final.taxa[1])
  orig.nrow<-nrow(dat_sp)
  nat_range <- iucn_anuran %>% filter(binomial %in% k) %>%
    st_union() %>% st_buffer(100) 
  dat_sp <- dat_sp %>% st_transform(crs.albers) %>% st_filter(nat_range)
  used.nrow <-nrow(dat_sp)
  nobstaxa<-data.frame(orig_nrow=orig.nrow,
             used_nrow=used.nrow,
             taxa=k,
             scientific_name=scientific_name)
  print(nobstaxa)
  counts_final<-bind_rows(counts_final, nobstaxa)
}
View(counts_final)
sum(counts_final$used_nrow)

