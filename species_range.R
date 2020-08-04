# Determining Anuran AOO within contiguous US - "species_range.R"
# edited on 2020Aug03 by TPD

# Install necessary libraries.
library(sf); library(tidyverse); library(scales)

# Set data paths
PATH_FocalSpecies_OccurrenceData <- "data/occ_data_used/"
PATH_tax_ref <- "data/AnuranTaxRef_20200708.csv"

# load the two spatial extents
US <- map_data("state")#converts state map to dataframe for mapping with ggplot2
usa <- maps::map("usa") #gets map for subsetting entries. Map of the lower 48 obtained through Package ‘maps’
#version 3.3.0 (Becker et al. 2018) in Program R (R Core Team 2016);
main <- map_data("usa", region=c("main"))#subsets map points to only mainland, no islands

#converts from dataframe to spatial polygon
df <- data.frame(main$long, main$lat) #dataframe of longitude and latitude 
usa_l48 <- st_as_sf(df, coords=c("main.long","main.lat"),
              crs=st_crs("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84")) %>%
  st_combine(.) %>%
  st_cast(., "POLYGON")

crs.geo <- st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84")
crs.albers<-st_crs("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs")
usa_l48_albs<-st_transform(usa_l48, crs.albers)

# Create empty AOO data table to be populated.
per_sp_range <- data.frame(scientific_name = character(), 
                           aoo_usal48 = numeric(),
                           aoo_na = numeric(),
                           stringsAsFactors = F)

# Create a list of occurrence data files for all focal species. 
tax_reference<-read.csv(PATH_tax_ref)
anuran.taxa<-unique(tax_reference$final.taxa)
head(anuran.taxa)

for(i in 1:length(anuran.taxa)){
  geodata<-read.csv(paste0(PATH_FocalSpecies_OccurrenceData, anuran.taxa[i], ".csv")) %>%
    dplyr::select(family, genus, species, Longitude, Latitude, year, source, final.taxa)
  
  # Convert from 'geodata' CSV to a sf object. 
  dat_sp <- st_as_sf(geodata, coords=c("Longitude","Latitude"), crs = crs.geo) #Save spatial dataframe with correct projection.
  
  # Set column 1 in AOOs to species name.
  per_sp_range[i,1] <- as.character(dat_sp$final.taxa[1])
  
  # Calculate area with 1km buffered circles inside the contiguous US 
  # we use the albers projection because it is projection
  dat_sp <- st_transform(dat_sp, crs.albers) #overwrites to the new projection (save on memory)
  
  # Create 1 km buffer and calculate the total area occupied per species.
  #   this number will be a result of the filtering completed at the occurrence data curation step
  dissolved_1km_buffer <- st_union(st_buffer(dat_sp, dist = 1))
  per_sp_range[i,]$aoo_na <- st_area(dissolved_1km_buffer)
  
  # Filter points that are within the contiguous US, create the buffer, & calculate area
  dat_sp_usa<-st_filter(dat_sp, st_as_sf(usa_l48_albs)) 
  d1km_buff_usa <- st_union(st_buffer(dat_sp_usa, dist = 1))
  per_sp_range[i,]$aoo_usal48 <- st_area(d1km_buff_usa)
}


# Calculate the percentage of area within the contiguous US
per_sp_range <- per_sp_range %>% select(-rowid) %>%
  mutate(per_inside_us=round((per_sp_range$aoo_usal48 / per_sp_range$aoo_na)*100, 2)) %>%
  arrange(per_inside_us) %>%
  rowid_to_column()
write.csv(per_sp_range, 'rcs_results/sp_range_within_US.csv')

#how many taxa would be removed if we only considered species with >20% of their range within the US
per_sp_range %>% arrange(per_inside_us) %>%
  filter(per_inside_us < 20) #7 species removed

ggplot()+
  stat_ecdf(data=per_sp_range, aes(x=per_inside_us))+
  scale_y_continuous(name="Proportion species")+
  scale_x_continuous(name="Percent Range within the Contiguous US")+
  theme_bw()+coord_flip()
ggsave('rcs_results/figures/sp_range_accum_curve.jpg',
       width=4, height=4)

per_sp_range_filt <- per_sp_range %>%
  filter(per_inside_us <50)
ggplot()+
  geom_line(data=per_sp_range_filt, aes(x=rowid, y=per_inside_us))+
  scale_x_continuous(name="Species",
                     breaks=per_sp_range_filt$rowid,
                     labels=per_sp_range_filt$scientific_name)+
  scale_y_continuous(name="Percentage Range within US")+
  theme_bw()+ 
  theme(axis.text.x = element_text(angle=30, hjust=1))
ggsave('rcs_results/figures/sp_range_accum_lowper.jpg',
       width=6, height = 4)

#### testing Smilisca baudinii and Dryophytes eximius ####
#Smilisca baudinii
sb_geodata<-read.csv(paste0(PATH_FocalSpecies_OccurrenceData, 
                            "Smilisca baudinii",
                            ".csv")) %>%
  dplyr::select(family, genus, species, Longitude, Latitude, year, source, final.taxa)
sb_dat_sp <- st_as_sf(sb_geodata, coords=c("Longitude","Latitude"), crs = crs.geo) #Save spatial dataframe with correct projection.
sb_dat_sp <- st_transform(sb_dat_sp, crs=crs.albers)
sb_dat_sp_usa<-st_filter(sb_dat_sp, st_as_sf(usa_l48_albs)) 

ggplot()+geom_sf(data=usa_l48_albs)+
  geom_sf(data=sb_dat_sp, color="blue")+
  geom_sf(data=sb_dat_sp_usa, color="red", alpha=0.3)
  
#Dryophytes eximius
de_geodata<-read.csv(paste0(PATH_FocalSpecies_OccurrenceData, 
                            "Dryophytes eximius",
                            ".csv")) %>%
  dplyr::select(family, genus, species, Longitude, Latitude, year, source, final.taxa)
de_dat_sp <- st_as_sf(de_geodata, coords=c("Longitude","Latitude"), crs = crs.geo) #Save spatial dataframe with correct projection.
de_dat_sp <- st_transform(de_dat_sp, crs=crs.albers)
de_dat_sp_usa<-st_filter(de_dat_sp, st_as_sf(usa_l48_albs)) 

ggplot()+geom_sf(data=usa_l48_albs)+
  geom_sf(data=de_dat_sp, color="blue")+
  geom_sf(data=de_dat_sp_usa, color="red", alpha=0.3)
