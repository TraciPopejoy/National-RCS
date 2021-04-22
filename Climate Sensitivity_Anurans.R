# National Fishes Vulnerability Assessment Project - "ARC_Climate Sensitivity_SCS.R"
# Revised by Sam Silknetter, 29June2020
# Revised by Traci DuBose, 2020Aug28

# This code generates standard deviations of climate variables.

# This code loads PRISM climate data and requires a HPC.
library(prism);library(tidyverse);library(raster);library(rgdal);library(sf); library(Hmisc)

# Set data paths for input data & load fixed functions ------

#where the all watershed shapefiles
PATH_HUC12 <- '/home/tracidubose/huc12_wgs84_sf_02010323.rds'
PATH_hydrosheds<-'/home/tracidubose/hydrobasins_wgs84_sf_20210323.rds'

#where species shapefiles reside
PATH_SHP_WS_L48 <- "/home/tracidubose/RCS_Results/20210330_sp_shp_files/huc12_sp/"
PATH_SHP_1km_L48 <- "/home/tracidubose/RCS_Results/20210330_sp_shp_files/L48_1km"
PATH_SHP_WS_NA <- "/home/tracidubose/RCS_Results/20210330_sp_shp_files/NA_hydrobasins/"
PATH_SHP_1km_NA <- "/home/tracidubose/RCS_Results/20210330_sp_shp_files/NA_1km"

#climate files
# PRISM files from Sam Silknetter (in press)
prism_ppt <- raster('/home/tracidubose/prism_files/ppt_mean_usa.tif')
prism_tmin<-raster("/home/tracidubose/prism_files/Tmin_mean_usa.tif")
prism_tmax <- raster("/home/tracidubose/prism_files/Tmax_mean_usa.tif") 
# loading all the WorldClim climate data files
NW_hemisphere<-extent(-180, -45, 0 ,90) #crop the raster to NW hemisphere
annual_ppt<-raster("/home/tracidubose/WorldClim_data/wc2.1_2.5m_bio_12.tif")
annual_ppt<-crop(annual_ppt, NW_hemisphere)
max_temp_sum<-raster("/home/tracidubose/WorldClim_data/wc2.1_2.5m_bio_10.tif")
max_temp_sum<-crop(max_temp_sum, NW_hemisphere)
min_temp_win<-raster("/home/tracidubose/WorldClim_data/wc2.1_2.5m_bio_11.tif")
min_temp_win<-crop(min_temp_win, NW_hemisphere)
crs.worldclim<-st_crs(annual_ppt)
#where the files end up:
PATH_ws_climate<-'/home/tracidubose/RCS_Results/20210330_watershed_climate_means/'

# Create a list of input shapefiles for loop.
FILES_watershed_NA <- list.files(path=PATH_SHP_WS_NA, pattern = "\\.shp$")
FILES_buffer_NA <- list.files(path=PATH_SHP_1km_NA, pattern = "\\.shp$")
FILES_watershed_L48 <- list.files(path=PATH_SHP_WS_L48, pattern = "\\.shp$")
FILES_buffer_L48 <- list.files(path=PATH_SHP_1km_L48, pattern = "\\.shp$")
anuran.taxa <- sub(".shp", "", FILES_watershed_L48)

# Create Mean and Standard Deviation functions to be used in script.
mean_1 <- function(x, na.rm = TRUE) {
  mean(x, na.rm=TRUE)}

sd_1 = function(x, ...){
  v = var(x, ...)
  l = sum(!is.na(x))
  sqrt(v*(l-1)/l)}

sd_2 = function(x, ...){# SD function to exclude non-zero values for binary for loop.
  v = var(x[!x == 0], ...)
  l = sum(!is.na(x))
  sqrt(v*(l-1)/l)}

# Make a binary dataframe of occupied watersheds per species ----
# Load Watershed shapefiles 
HUC12_All <- readRDS(PATH_HUC12)
st_crs(HUC12_All)

hucIDlong<-NULL
for(i in 1:length(FILES_watershed_L48)){
  #read in the occupied watersheds
  shp<-st_read(PATH_SHP_WS_L48 , sub('.shp', '', FILES_watershed_L48[i]),
               quiet = T)
  #rearrange data to a long list of occupied watersheds
  hucIDlong.one<-data.frame(taxa=rep(sub(".shp", "", FILES_watershed_L48[i]), 
                                     length(shp$wtrshd_)),
                            watershed_ID=shp$wtrshd_) #%>%
    #mutate_if(is.factor, as.character)
  #add it all together into one long list
  hucIDlong<-rbind(hucIDlong, hucIDlong.one)
}

# reducing huc spatial file for easier processing
huc12_red<-HUC12_All %>%
  filter(huc12 %in% unique(hucIDlong$watershed_ID))

#convert this to a binary matrix for easy multiplication later
hucIDbin <- hucIDlong %>% 
  mutate(watershed_ID = as.character(paste(watershed_ID)),
         presence=1) %>%
  pivot_wider(names_from='taxa', values_from = 'presence', values_fill=list(presence = 0))
write.csv(hucIDbin, '/home/tracidubose/RCS_Results/HUC12_binary_matrix_20210330.csv', row.names = F)

hydroIDlong<-NULL
for(i in 1:length(FILES_watershed_NA)){
  shp<-st_read(PATH_SHP_WS_NA , sub('.shp', '', FILES_watershed_NA[i]),
               quiet = T)
  hydroIDlong.one<-data.frame(taxa=rep(sub(".shp", "", FILES_watershed_NA[i]), 
                                     length(shp$wtrshd_)),
                            watershed_ID=shp$wtrshd_) #%>%
  #mutate_if(is.factor, as.character)
  hydroIDlong<-rbind(hydroIDlong, hydroIDlong.one)
}

# reducing huc spatial file for easier processing
Hydro_all <- readRDS(PATH_hydrosheds)
st_crs(Hydro_all)
hydro_red<-Hydro_all %>%
  filter(HYBAS_ID %in% unique(hydroIDlong$watershed_ID))
#ggplot()+geom_sf(data=hydro_red)

hydroIDbin <- hydroIDlong %>% 
  mutate(presence=1) %>%
  pivot_wider(names_from='taxa', values_from = 'presence', 
              values_fill=list(presence = 0))
write.csv(hydroIDbin, '/home/tracidubose/RCS_Results/Hydro12_binary_matrix_20210330.csv', row.names = F)


#### Watershed Climate Mean Calculations ####
# I use the extract_mean_env function to pull the weighted mean climate values from each watershed.
# this can take a while (>24 hrs on ARC). As such, I run them once and write out the results.
# Function to extract and quantify an environmental variable from a raster by watershed 
extract_mean_env <- function(env.raster, 
                             hucs, 
                             ID.col,
                             hucbinary,
                             write.file=T, 
                             path='/home/tracidubose/prism_files/'){
  
  hucs <- hucs %>% rename(watershed.name=all_of(ID.col)) #make this easier for later
  
  # extract the values of interest from the raster for each polygon
  env.values<-raster::extract(env.raster, hucs, 
                              fun = mean, na.rm = T, weights = TRUE, 
                              small = TRUE, method = 'simple', df = TRUE)
  env.values$watershed.name <- hucs$watershed.name # add id column back to resulting dataframe
  env.values <- env.values[, c(3, 2)] #reorder resulting matrix and lose ID column
  if(write.file==T){write.csv(env.values, paste0(path, ID.col, "_extracted_", 
                                                  deparse(substitute(env.raster)),
                                                  '.csv'))} #writing out a csv just in case
  return(env.values)
}

# running the function with either the hydrobasins + worldclim or huc12 + prism
# ONLY NEED TO RUN THIS ONCE! THE HUCS TAKE A LONG TIME
WS_NA_ppt<-extract_mean_env(env.raster=annual_ppt, 
                            hucs=hydro_red, 
                            ID.col="HYBAS_ID",
                            write.file=T, 
                            path=PATH_ws_climate)
WS_NA_tmax<-extract_mean_env(env.raster=max_temp_sum, 
                            hucs=hydro_red, 
                            ID.col="HYBAS_ID",
                            write.file=T, 
                            path=PATH_ws_climate)
WS_NA_tmin<-extract_mean_env(env.raster=min_temp_win, 
                            hucs=hydro_red, 
                            ID.col="HYBAS_ID",
                            write.file=T, 
                            path=PATH_ws_climate)

WS_L48_ppt<-extract_mean_env(env.raster=prism_ppt, 
                             hucs=huc12_red, 
                            ID.col="huc12",
                            write.file=T, 
                            path=PATH_ws_climate)
WS_L48_tmax<-extract_mean_env(env.raster=prism_tmax, 
                             hucs=huc12_red, 
                             ID.col="huc12",
                             write.file=T, 
                             path=PATH_ws_climate)
WS_L48_tmin<-extract_mean_env(env.raster=prism_tmin, 
                              hucs=huc12_red, 
                              ID.col="huc12",
                              write.file=T, 
                              path=PATH_ws_climate)

# pull in the extracted climate means
for(k in list.files(PATH_ws_climate, pattern='extract')){
  assign(sub('extracted_', '', sub('.csv', '', k)), 
         read.csv(paste0(PATH_ws_climate, k),
                  row.names = 1, stringsAsFactors = F) %>%
           mutate(watershed.name=ifelse(nchar(watershed.name)==11, 
                                      paste0('0', watershed.name), watershed.name)))
} 

#check the watershed IDs still match (if character ->, often drop the leading 0)
which(nchar(huc12_prism_ppt$watershed.name)==11)
which(!(hydroIDbin$watershed_ID %in% HYBAS_ID_annual_ppt$watershed.name))

# use matrix algebra to build a matrix to calculate standard of deviation
# the binary matrix should contain 0 where sp is not present and mean environmental values where sp is present
# need to turn that into a long dataframe that excludes 0s 
hucIDbin_ordered<-hucIDbin[order(hucIDbin$watershed_ID),] %>%
  dplyr::select(-watershed_ID)
hydroIDbin_ord<-hydroIDbin[order(hydroIDbin$watershed_ID),] %>%
  dplyr::select(-watershed_ID)
#list of huc12 climate dataframes
climate_means<-sub('extracted_', '', sub('.csv', '', list.files(PATH_ws_climate, pattern='extract'))) 
clim_var<-c('mean','sd')
for(i in climate_means){
  climate_df<-get(i)
  if(substr(i,1,3)=="huc"){
    climate_df<-climate_df %>% 
      filter(watershed.name %in% hucIDbin$watershed_ID)
    print(nrow(hucIDbin_ordered) == nrow(climate_df)) #check to make sure number of rows match
    climate_binmat<-climate_df[order(climate_df$watershed.name),2] * hucIDbin_ordered
  }
  if(substr(i,1,3)=='HYB'){
    print(nrow(hydroIDbin_ord) == nrow(climate_df))
    climate_binmat<-climate_df[order(climate_df$watershed.name),2] * hydroIDbin_ord
  }
  assign(paste0(i, '_taxa_sum'),
    climate_binmat %>%
    as_tibble() %>%
    bind_cols(watershed.name=climate_df[order(climate_df$watershed.name),1]) %>%
    pivot_longer(-watershed.name, names_to='species') %>%
    filter(value!=0) %>%
    mutate(watershed.name=as.character(watershed.name),
           value_origin=substr(i,1,5),
           value_type=case_when(grepl('ppt', names(climate_df)[2])| grepl('bio_12', names(climate_df)[2]) ~ 'annual ppt',
                                grepl('Tmax', names(climate_df)[2])| grepl('bio_10', names(climate_df)[2]) ~ 'tmax',
                                grepl('Tmin', names(climate_df)[2])| grepl('bio_11', names(climate_df)[2]) ~ 'tmin')))
}

head(huc12_prism_ppt_taxa_sum)
head(HYBAS_ID_annual_ppt_taxa_sum)

# Mean and Standard Deviation Calculation
taxa_climate_values_raw<-bind_rows(mget(paste0(climate_means, '_taxa_sum')))
write.csv(taxa_climate_values_raw, paste0('/home/tracidubose/RCS_Results/taxa_ws_climate_values',
                                   format(Sys.Date(), "%Y%m%d"),'.csv'))
#supplementary information
taxa_climate_values_raw %>%
  group_by(species, value_origin) %>% tally() %>%
  mutate(n=n/3) %>%
  pivot_wider(names_from=value_origin, values_from=n)

taxa_climate_values_raw %>%
  group_by(species, value_type) %>%
  dplyr::summarize(mean=mean(value), sd=sd(value))

taxa_climate_values_raw %>%
  group_by(value_type, value_origin) %>%
  ggplot()+
  geom_density(aes(x=value, fill=value_origin), alpha=.5)+
  #geom_vline(aes(xintercept=mean(value), group=value_type))+
  facet_wrap(~value_type, scales='free')+
  theme_bw() + theme(#axis.text.x = element_text(angle=20),
                     axis.title.x=element_text(size=-1),
                     legend.position = 'bottom')
ggsave('/home/tracidubose/RCS_Results/sup_fig/ws_climate_var_dist_20210402.jpg',
       width=6, height=3)

#### Buffered Pts Climate Mean Calculations ####
#Extract environmental data from each species 1km buffer AOO. 
extract_mean_env_buff <- function(env.raster, 
                                  variable,
                                  spp, 
                                  native_range=F,
                                  buff_file_folder,
                                  write.f=T,
                                  output_path){
  
  
  if(length(list.files(buff_file_folder, pattern=spp))!=0){
    spp_BUF <- read_sf(buff_file_folder, spp) #load the buffered shape file
    if(native_range==T){
      spp_BUF<-st_intersection(iucn_anuran[iucn_anuran$binomial==spp,], spp_BUF) %>%
        st_union() %>% st_as_sf()
    }
    spp_BUF <-  st_transform(spp_BUF, st_crs(env.raster)) #transform to match raster
  # extract the values of interest from the raster for each polygon
  env.values<-raster::extract(env.raster, spp_BUF, 
                              na.rm = T, weights = TRUE, 
                              small = TRUE, method = 'simple', df = TRUE)
  names(env.values)[2]<-variable
  source<-substr(deparse(substitute(env.raster)),1,5)
  if(write.f==T){write.csv(env.values %>% mutate(species=spp,
                                                 source=source),
                           paste0(output_path, spp, '_', deparse(substitute(env.raster)),'.csv'))}
  print(paste(spp, deparse(substitute(env.raster))))
  }else{
    print(paste(spp, 'has no area outside of the contiguous US'))
  }
}

st_transform(read_sf(PATH_SHP_1km_NA, 'Anaxyrus canorus'), st_crs(prism_ppt))

crs.albers <- st_crs("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs")#2.12kb
iucn_anuran<-read_sf('/home/tracidubose/RCS_Anuran_input/ANURA/', 'ANURA') %>% 
  st_transform(., crs.albers)
st_crs(iucn_anuran)
#note: when native = T, will get a lot of warnings for two taxa with native ranges not reported at IUCN:
#Acris blanchardi & Lithobates kauffauldi

PATH_buff_climate_val<-'/home/tracidubose/RCS_Results/20210330_watershed_climate_means/buffer climates/'
taxa.toNOTdo<-gsub('\\_.*','', list.files(PATH_buff_climate_val)) %>%
  as_tibble() %>%
  count(value) 

taxa.toNOTdo %>% filter(!(n %in% c(6,3)))
taxa.todo<-anuran.taxa[!(anuran.taxa %in% taxa.toNOTdo$value)]
taxa.todo

grep('woodhou', list.files(PATH_buff_climate_val), value=T)

#NOTE: Anaxyrus baxteri and Anaxyrus canorus messed up - outside native range?
crs.geo <- st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84")
read.csv('/home/tracidubose/RCS_Anuran_input/Anuran Occurrence Data/Anaxyrus canorus_20210330.csv') %>% 
  st_as_sf(coords=c("Longitude","Latitude"), crs=crs.geo) %>%
  st_transform(st_crs(iucn_anuran))%>% 
  ggplot()+
  geom_sf(data=iucn_anuran %>% filter(binomial=="Anaxyrus canorus"))+
  geom_sf()
#yep, need to go back and fix by hand? or just ignore?
list.files(PATH_SHP_1km_NA, pattern='canor')


for(u in taxa.todo){
  extract_mean_env_buff(env.raster = annual_ppt, variable='ppt', spp=u, native_range = F,
                        buff_file_folder = PATH_SHP_1km_NA, output_path = PATH_buff_climate_val)
  extract_mean_env_buff(env.raster = prism_ppt, variable='ppt', spp=u, native_range = F,
                        buff_file_folder = PATH_SHP_1km_L48, output_path = PATH_buff_climate_val)
  
  extract_mean_env_buff(env.raster = min_temp_win, variable='tmin', spp=u, native_range = F,
                        buff_file_folder = PATH_SHP_1km_NA, output_path = PATH_buff_climate_val)
  extract_mean_env_buff(env.raster = prism_tmin, variable='tmin', spp=u, native_range = F,
                        buff_file_folder = PATH_SHP_1km_L48, output_path = PATH_buff_climate_val)
  
  extract_mean_env_buff(env.raster = max_temp_sum, variable='tmax', spp=u, native_range = F,
                        buff_file_folder = PATH_SHP_1km_NA, output_path = PATH_buff_climate_val)
  extract_mean_env_buff(env.raster = prism_tmax, variable='tmax', spp=u, native_range = F,
                        buff_file_folder = PATH_SHP_1km_L48, output_path = PATH_buff_climate_val)
  
}


clim_values<-NULL
for(b in list.files(PATH_buff_climate_val,
                    full.name=T)){
  ms<-read.csv(b)
  clim_values<-bind_rows(clim_values, ms)
}
which(!(unique(clim_values$species) %in% anuran.taxa)) #wohoo! have all taxa done
clim_values %>% dplyr::select(-X,-ID) %>%
  as_tibble() %>% 
  pivot_longer(cols = c('ppt','tmin', 'tmax')) %>%
  filter(!is.na(value)) %>%
  group_by(species, source, name) %>%
  dplyr::summarize(sumweight=sum(weight)) #some of them have weights lower than 1 (NA values)
  
(aoos<-read.csv('/home/tracidubose/RCS_Results/AOOs_raw_20210331.csv') %>%
  filter(!(area.sqkm==0 & scientific_name %in% c('Acris blanchardi','Lithobates kauffeldi')),
         grepl('1km', area.type)) %>%
  group_by(scientific_name) %>%
  mutate(prop.area=area.sqkm/sum(area.sqkm)) %>%
  dplyr::select(-X))

buf_climate_sum<-clim_values %>% 
  dplyr::select(-X,-ID) %>%
  as_tibble() %>% 
  pivot_longer(cols = c('ppt','tmin', 'tmax')) %>%
  filter(!is.na(value)) %>%
  group_by(species, source, name) %>%
  mutate(area.type=case_when(source=='prism'~'usa_l48_alb_1km',
                             source!='prism'~'NOTusa_l48_alb_1km')) %>%
  left_join(aoos, by=c('species'='scientific_name', 'area.type')) %>%
  mutate(new.weight=weight*prop.area) %>%
  group_by(species, name, source) %>%
  dplyr::summarize(mean_value=wtd.mean(value, new.weight, na.rm = TRUE, normwt = TRUE),
                   sd_value=sqrt(wtd.var(value, new.weight, na.rm = TRUE, normwt = TRUE)),
                   sumweight=sum(new.weight))

write.csv(buf_climate_sum, 
          paste0('/home/tracidubose/RCS_Results/20210330_watershed_climate_means/buff_climate_sums_',
                 format(Sys.Date(), '%Y%m%d'),'.csv'))

range.wide.weighted.CS <- buf_climate_sum %>% 
  group_by(species, name) %>%
  dplyr::summarise(mean_value=wtd.mean(mean_value, sumweight, na.rm = TRUE, normwt = TRUE),
                   sd_value=wtd.mean(sd_value, sumweight, na.rm = TRUE, normwt = TRUE))
write.csv(range.wide.weighted.CS, 
          paste0('/home/tracidubose/RCS_Results/20210330_watershed_climate_means/range.wide.weighted.CS_',
                 format(Sys.Date(), '%Y%m%d'),'.csv'))
