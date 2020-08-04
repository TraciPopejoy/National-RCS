# Area of Occupancy Code for the National Anuran RCS Assessment - "Area of Occupancy Code_Anurans.R"
# Revised by Sam Silknetter, 03April2020
# edited on 2020July18 by TPD

# Install necessary libraries.
library(sf); library(tidyverse); library(scales)

# Set data paths.
PATH_HUC12 <- '/home/tracidubose/huc12_wgs84_sf.rds'
PATH_USAL48<- '/home/tracidubose/usal48_nad83_sf.rds'
PATH_FocalSpecies_OccurrenceData <- "/home/tracidubose/RCS_Anuran_input/Anuran Occurrence Data/"
PATH_SHP_HUC12 <- "/home/tracidubose/rcs_results/huc12_sp/"
PATH_SHP_dis1km <- "/home/tracidubose/rcs_results/dis1km/"

# Set spatial coordinate reference system (CRS).
crs.geo <- st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84")
crs.albers <- st_crs("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs")#2.12kb

# Read in fixed data.
huc12<-readRDS(PATH_HUC12)
st_crs(huc12)==crs.geo #checking projection is correct
usa_l48_albs<-readRDS(PATH_USAL48)
st_crs(usa_l48_albs)==crs.albers

# Create empty AOO data table to be populated.
AOOs <- data.frame(scientific_name = character(), dissolved_buffer_1km = numeric(), 
        huc12_area = numeric(), 
        mean_Rank = numeric(), range = numeric(), sd_Rank = numeric(), StandardizedMean = numeric(), 
        cv = numeric(), logBuffer1km = numeric(),
        logHUC12 = numeric(), DFdif=numeric(),
        stringsAsFactors = F)

# Create a list of occurrence data files for all focal species. 
tax_reference<-read.csv("/home/tracidubose/RCS_Anuran_input/AnuranTaxRef_20200708.csv")
anuran.taxa<-unique(tax_reference$final.taxa)
head(anuran.taxa)

for(i in 1:length(anuran.taxa)){
  #load the data into geodata [!not split by time yet]
  geodata<-read.csv(paste0(PATH_FocalSpecies_OccurrenceData, anuran.taxa[i], ".csv")) %>%
      dplyr::select(family, genus, species, Longitude, Latitude, year, source, final.taxa)
    
  # Convert from 'geodata' CSV to a sf object. 
  dat_sp <- st_as_sf(geodata, coords=c("Longitude","Latitude"), crs = crs.geo) #Save spatial dataframe with correct projection.
    
  # Set column 1 in AOOs to species name.
  AOOs[i,1] <- as.character(dat_sp$final.taxa[1])

  # Calculate area (km^2) of occupied watersheds per species. 
  #! need to make sure these name columns match whichever spatial dataframe I'm using (line 62 & line 65)
  dat_sp <- st_join(dat_sp, huc12) # Store the HUC12 name as an attribute of the fish data.
  #note that this can add columns
  AOOs[i,]$DFdif<-nrow(geodata)-nrow(dat_sp)
  Number_of_HUC12s <- data.frame(huc12=unique(dat_sp$huc12))
  huc12_sp <- huc12[huc12$huc12 %in% Number_of_HUC12s$huc12,]  %>%
    select(-noncontributingareaacres, -noncontributingareasqkm)    # Extract unique HUC12 data for the species.
  AOOs[i,]$huc12_area <- sum(dplyr::select(as.data.frame(huc12_sp), contains('sqkm')), na.rm=T)  # Sum the total area (square kilometer) for each unique HUC12.
  
  # Generate 1km point buffers. You must use a projected CRS for function gBuffer(). For CRS, we use: USA Contiguous albers equal area. 
  dat_sp <- st_transform(dat_sp, crs.albers) #overwrites to the new projection (save on memory)
    
  # Create 1 km buffer and calculate the total area occupied per species.
  dat_sp_us<-st_filter(dat_sp, st_as_sf(usa_l48_albs)) 
  dissolved_1km_buffer <- st_union(st_buffer(dat_sp_us, dist = 1))
  AOOs[i,]$dissolved_buffer_1km <- st_area(dissolved_1km_buffer)
    
  # Save the spatial dataframes as ESRI Shapefiles. 
  species <- geodata$final.taxa[1] # Identify the species name for CSV output. 
  st_write(huc12_sp, dsn = paste0(PATH_SHP_HUC12, species, "_HUC12.shp"))
  st_write(dissolved_1km_buffer, dsn = paste0(PATH_SHP_dis1km, species, "_1km.shp"))
}


# Calculate AOO rank for each range metric and focal species.
AOOs$rank1km <- rank(AOOs$dissolved_buffer_1km)
AOOs$rank_huc12 <- rank(AOOs$huc12_area)

for(i in 1:length(anuran.taxa)){
  # Calculate and store mean of the 2 AOO ranks in the mean column.
  AOOs[i,]$mean_Rank <- mean(c(AOOs[i,]$rank1km, AOOs[i,]$rank_huc12), na.rm=TRUE)

  # Calculate and store the range of the 2 AOO ranks in the range column.
  AOOs[i,]$range <- abs((AOOs[i,]$rank1km - AOOs[i,]$rank_huc12 )) #$rank1km, AOOs[i,]$rank_huc12), na.rm=TRUE)

  # Calculate and store standard deviation of the 2 AOO ranks in the sd column.
  AOOs[i,]$sd_Rank <- sd(c(AOOs[i,]$rank1km, AOOs[i,]$rank_huc12), na.rm=TRUE)

  # Scale mean rank between 0 and 1
  AOOs$StandardizedMean <- rescale(AOOs$mean_Rank, to = c(0,1))

  # Calculate and store coefficient of variation (cv) for rankings.
  AOOs$cv <- (AOOs$sd_Rank/AOOs$mean_Rank)*100

  # Calculate and store logarithms (base 10) of AOOs for both range metrics.
  AOOs$logBuffer1km <- log10(AOOs$dissolved_buffer_1km)
  AOOs$logHUC12 <- log10(AOOs$huc12_area)
}

# Write CSV with date in "_01Jan2020" format.
write.csv(AOOs, file = paste0("rcs_results/AOO HUC12 Output_",
                              format(Sys.Date(), '%Y%m%d'),'.csv'))
