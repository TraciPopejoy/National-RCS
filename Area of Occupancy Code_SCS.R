# Area of Occupancy Code for the National Fishes Vulnerability Assessment - "Area of Occupancy Code_SCS.R"
# Revised by Sam Silknetter, 03April2020
# edited on 2020July01 by TPD

### NOTES: need to download HUC12 shapefile - look for original source to match Sam&Jen&Abby
###       currently have all HUC12 lines commented out
###       if we want to look at rcs through time, might want to:
                # 1. pull in data
                # 2. apply time breaks
                # 3. build AOOs & CS score based on the smaller sample (filter function)
                # 4. AOO & CS table will both species & period identifier

# Install necessary libraries.
library(raster);library(rgdal);library(sp);library(rgeos);library(scales)

# Set data paths.
#PATH_HUC12 <- "G:/Shared drives/NFVAP - National RCS/R Code/National-RCS/HUC12"
  # Use only the subset of candidate spp. that met filtering requirements or are 'exception species'. These are the 'focal species'.  
PATH_FocalSpecies_OccurrenceData <- "occ_data"
#PATH_SHP_HUC12 <- "G:/Shared drives/NFVAP - National RCS/R Code/National-RCS/Shapefiles/HUC-12/"
PATH_SHP_1km <- "rcs_results/1kmbuf/"
PATH_SHP_dis1km <- "rcs_results/dis1km/"

# Set spatial coordinate reference system (CRS).
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
crs.albers <- CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs")

# Read in fixed data.
#huc12 <- readOGR(dsn = PATH_HUC12, layer = "HUC12")
#proj4string(huc12) <- CRS("+init=epsg:4269 +proj=longlat +ellps=GRS80 +datum=NAD83")
#huc12_wgs84 <- spTransform(huc12, crs.geo)

# Create empty AOO data table to be populated.
AOOs <- data.frame(scientific_name = character(), buffer_1km = numeric(), dissolved_buffer_1km = numeric(), 
        huc12_area = numeric(), mean_Rank = numeric(), range = numeric(), sd_Rank = numeric(), StandardizedMean = numeric(), 
        cv = numeric(), logBuffer1km = numeric(), logHUC12 = numeric(), stringsAsFactors = F)

# Create a list of occurrence data files for all focal species. 
#FILES <- list.files(path=PATH_FocalSpecies_OccurrenceData)
source('LoadingOccDataframes.R') #load the occurrence datasets
anuran_taxa<-unique(Anurans_df$species)

# Loop to cycle through all species
orig.time<-paste(Sys.time())
#for(i in 1:length(unique(anuran_taxa))){ 
  
for(i in 1:5){
  #problem: some species have multiple names within their datasheet - this causes problems!
  
  
  #load the data into geodata [!not split by time yet]
  geodata <- Anurans_df %>% 
    filter(species==unique(anuran_taxa)[i]) %>% # only pull out one taxa (matches i)
    dplyr::select(family, genus, species, Longitude, Latitude, year, source, #important columns
           gbifID, ID, UUID)#keeping these columns just in case

  # Convert from 'geodata' CSV to an OGR object. 
  dat_csv <- as.data.frame(geodata)
  dat_coord <- dat_csv[c("Longitude", "Latitude")] # Create data frame of longitude and latitude for all entries.
  dat_sp <- SpatialPointsDataFrame(dat_coord, dat_csv, proj4string = crs.geo) # Save spatial dataframe with correct projection.
  
  # Set column 1 in AOOs to species name.
  AOOs[i,1] <- as.character(dat_sp$species[1])  
  
  # Calculate area (km^2) of occupied watersheds per species. 
  #dat_sp$HUC12 <- over(dat_sp, huc12_wgs84)$HUC12 # Store the HUC12 name as an attribute of the fish data.
  #Number_of_HUC12s <- data.frame(unique(dat_sp$HUC12))
  #huc12_sp <- huc12[huc12$HUC12 %in% Number_of_HUC12s$unique.dat_sp.HUC12.,]  # Extract unique HUC12 data for the species.
  #total_area_huc12 <- sum(huc12_sp$AreaSqKm) # Sum the total area (square kilometer) for each unique HUC12.
  #AOOs[i,]$huc12_area <- total_area_huc12
  
  # Generate 1km point buffers. You must use a projected CRS for function gBuffer(). For CRS, we use: USA Contiguous albers equal area. 
  geodata_Albers <- spTransform(dat_sp, crs.albers)
  
  # Create 1 km buffer and calculate the total area occupied per species.
  sp_buffer_1km <- gBuffer(geodata_Albers, width = 1, byid= TRUE)
  AOOs[i,]$buffer_1km <- gArea(sp_buffer_1km, byid = FALSE) #to be deleted, but keeping for now as a check. 
  dissolved_1km_buffer <- gUnaryUnion(sp_buffer_1km, id=NULL)
  AOOs[i,]$dissolved_buffer_1km <- gArea(dissolved_1km_buffer, byid = FALSE)
  
  # Save the spatial dataframes as ESRI Shapefiles. 
  species <- geodata$species[1] # Identify the species name for CSV output. 
  #raster::shapefile(huc12_sp, file = paste0(PATH_SHP_HUC12, species, "_HUC12.shp"), overwrite=TRUE)
  raster::shapefile(sp_buffer_1km, file = paste0(PATH_SHP_1km, species, "_1km.shp"), overwrite=TRUE)
  raster::shapefile(dissolved_1km_buffer, file = paste0(PATH_SHP_dis1km, species, "_1km.shp"), overwrite=TRUE)
}


# Calculate AOO rank for each range metric and focal species.
AOOs$rank1km <- rank(AOOs$dissolved_buffer_1km)
#AOOs$rank_huc12 <- rank(AOOs$huc12_area)

for(i in 1:30){
  # Calculate and store mean of the 2 AOO ranks in the mean column.
  AOOs[i,]$mean_Rank <- mean(c(AOOs[i,]$rank1km, AOOs[i,]$rank_huc12), na.rm=TRUE)

  # Calculate and store the range of the 2 AOO ranks in the range column.
  AOOs[i,]$range <- abs((AOOs[i,12] - AOOs[i,13])) #$rank1km, AOOs[i,]$rank_huc12), na.rm=TRUE)

  # Calculate and store standard deviation of the 2 AOO ranks in the sd column.
  AOOs[i,]$sd_Rank <- sd(c(AOOs[i,]$rank1km, AOOs[i,]$rank_huc12), na.rm=TRUE)

  # Scale mean rank between 0 and 1
  AOOs$StandardizedMean <- rescale(AOOs$mean_Rank, to =c(0,1))

  # Calculate and store coefficient of variation (cv) for rankings.
  AOOs$cv <- (AOOs$sd_Rank/AOOs$mean_Rank)*100

  # Calculate and store logarithms (base 10) of AOOs for both range metrics.
  AOOs$logBuffer1km <- log10(AOOs$dissolved_buffer_1km)
  #AOOs$logHUC12 <- log10(AOOs$huc12_area)
  print(AOOs[i,]$scientific_name)
}

# Write CSV with date in "_01Jan2020" format.
write.csv(AOOs, file = paste0("rcs_results/AOO Output_",
                              format(Sys.Date(), '%Y%m%d'),'.csv'))
