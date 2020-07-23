# National Fishes Vulnerability Assessment Project - "ARC_Climate Sensitivity_SCS.R"
# Revised by Sam Silknetter, 29June2020
# Revised by Traci DuBose, 19July2020

# This code generates standard deviations of climate variables.

# This code loads PRISM climate data and requires a HPC.
library(prism);library(tidyverse);library(raster);library(rgdal);library(sf)

# Set data paths for input data & load fixed functions ------
PATH_SHP_HUC12_AOO <- "/home/tracidubose/rcs_results/huc12_sp"
PATH_SHP_1km <- "/home/tracidubose/rcs_results/dis1km"
PATH_HUC12 <- "/home/tracidubose/huc12_nad83_sf.rds"
PATH_PRISM_ppt <- "/home/tracidubose/prism_files/ppt/"
PATH_PRISM_tmax <- "/home/tracidubose/prism_files/Tmax/"
PATH_PRISM_tmin <- "/home/tracidubose/prism_files/Tmin/"
PATH_PRISM_tmaxAug <- "/home/tracidubose/prism_files/TmaxAug/"
PATH_PRISM_tminJan <- "/home/tracidubose/prism_files/TminJan/"

# Create a list of input shapefiles for loop.
FILES_watershed <- list.files(path=PATH_SHP_HUC12_AOO, pattern = "\\.shp$")
FILES_buffer <- list.files(path=PATH_SHP_1km, pattern = "\\.shp$")

# Load Watershed shapefiles (as .rds) from path.
HUC12_All <- readRDS(PATH_HUC12)
st_crs(HUC12_All)
#loading in the base version - crs is NAD83/ epsf:4269

# Make a binary dataframe of occupied watersheds per species 
#!skip if already done once
#hucIDlong<-NULL
#for(i in 1:length(FILES_watershed)){
#  shp<-st_read(paste0(PATH_SHP_HUC12_AOO, "/", FILES_watershed[i]),
#               quiet = T)
#  hucIDlong.one<-data.frame(taxa=rep(sub("_HUC12.shp", "", FILES_watershed[i]), 
#                                     length(shp$huc12)),
#                            HUC_ID=shp$huc12) %>%
#    mutate_if(is.factor, as.character)
#  hucIDlong<-rbind(hucIDlong, hucIDlong.one)
#}

#check it got all taxa, also count number of hucs for each taxa
#hucIDlong %>% group_by(taxa) %>% tally() %>% arrange(n)

#length(unique(hucIDlong$HUC_ID))
#hucIDbin<-hucIDlong %>% mutate(present=1,
#                               HUC_ID=as.character(HUC_ID)) %>%
#  pivot_wider(names_from=taxa, values_from = present, values_fill=list(present=0)) %>%
#  arrange(HUC_ID) 

#which(is.na(hucIDbin$HUC_ID), arr.ind = T) #check for NA Huc IDs

#write.csv(hucIDbin, file = paste0("/home/tracidubose/rcs_results/OccupiedHuc12_",
#                                  format(Sys.Date(), "%Y%m%d"), ".csv"), row.names=FALSE)

# Import HUC12 x Anuran species presence matrix (binary matrix where 1= species occurrence)
# note! this doesn't preserve huc12 as a factor or character
hucIDbin<-read.csv("/home/tracidubose/rcs_results/OccupiedHuc12_20200721.csv",
                   colClasses = c('HUC_ID'="character"))
which(nchar(hucIDbin$HUC_ID)!=12)

#this data frame is used to make the groupings for running the extract_mean_env function piecemeal
group_hucs<-hucIDbin %>% dplyr::select(HUC_ID) %>%
  mutate(HUC2=substr(HUC_ID,1,4)) %>% 
  group_by(HUC2) %>%
  tally() %>% arrange(desc(n))
nrow(group_hucs)

# reducing huc spatial file for easier processing
huc12_red<-HUC12_All[match(hucIDbin$HUC_ID, HUC12_All$huc12),]
#pdf('/home/tracidubose/AnuranHucCheck.pdf')
#plot(huc12_red)
#dev.off()

# checking to see if order of hucIDbin is order of huc12_red
#head(hucIDbin$HUC_ID, n=7)
#head(as.character(huc12_red$huc12), n=7)
#hucIDbin$HUC_ID == huc12_red$huc12
#which(hucIDbin$HUC_ID != huc12_red$huc12) #they all match!

# Create Mean and Standard Deviation functions to be used in script.
mean_1 <- function(x, na.rm = TRUE) {
  mean(x, na.rm=TRUE)}

sd_1 = function(x, ...){
  v = var(x, ...)
  l = sum(!is.na(x))
  sqrt(v*(l-1)/l)}

  # SD function to exclude non-zero values for binary for loop.
sd_2 = function(x, ...){
  v = var(x[!x == 0], ...)
  l = sum(!is.na(x))
  sqrt(v*(l-1)/l)}

# Function to extract and quantify an environmental variable from a raster by hucs
extract_mean_env <- function(env.raster, hucs, write=T, 
                             path='/home/tracidubose/prism_files/', hucbinary){
  
  # extract the values of interest from the raster for each polygon
  env.values<-raster::extract(env.raster, hucs, fun = mean,na.rm = T, weights = TRUE, 
                              small = TRUE, method = 'simple', df = TRUE)
  env.values$HUC_ID <- hucs$huc12 # add huc_id column back to resulting dataframe
  # don't need below code if pulling from a saved raster
  #names(env.values)[names(env.values) == "layer"] <- paste0("mean_" , deparse(substitute(env.raster))) #rename layer to appropriate evironmental variable
  env.values <- env.values[, c(3, 2)] #reorder resulting matrix and lose ID column
  if(write==T){write.csv(env.values, paste0(path, "extract_output_", 
                                            deparse(substitute(env.raster)),
                                            '.csv'))} #writing out a csv just in case
  
  # check that each matrix is ordered the same (hucs & rownumbers aligned)
  test_order <- as.character(env.values$HUC_ID) == hucbinary$HUC_ID
  if(length(test_order) > sum(test_order)){ print("huc order not aligned")
    hucbinary<-hucbinary[match(env.values$HUC_ID, hucbinary$huc12),]}else{print("huc aligned")}
  
  # Use matrix algebra to create a dataframe of environmental variable per species.
  output <- hucbinary[,-1] * c(env.values[,2])
  row.names(output)<-hucbinary$HUC_ID
  return(output)
}

#### PRECIPITATION CALCS ####

# Set path and load 'raster stack' to variable-specific prism data. 
#options(prism.path = PATH_PRISM_ppt) 
#PRISM_ppt <- prism_stack(ls_prism_data()) # Note: Ensure prism.path is set to correct climate variable

# Mean and Standard Deviation Calculation
#mean_ppt <- raster::calc(PRISM_ppt, fun = mean_1, na.rm = TRUE)
#sd_ppt <- raster::calc(PRISM_ppt, fun = sd_1, na.rm=TRUE)
mean_ppt<-raster('/home/tracidubose/prism_files/ppt_mean_usa.tif')
sd_ppt<-raster('/home/tracidubose/prism_files/ppt_sd_usa.tif')

# I can't do it all at one go - it takes over 5 hrs for the extraction.
# So, I am doing it piecemeal with the following for loops

ppt.huc12.matrix<-NULL
for(u in group_hucs$HUC2[151:211]){
  hucsub<-hucIDbin %>%
    filter(substr(HUC_ID,1,4)==u)
  
  huc.sf.red<-huc12_red[match(hucsub$HUC_ID, huc12_red$huc12),]
  which(hucsub$HUC_ID != huc.sf.red$huc12)
  
  ppt.huc12.matrix.sub<-extract_mean_env(env.raster=mean_ppt, 
                                     hucs=huc.sf.red, 
                                     hucbinary=hucsub)
  ppt.huc12.matrix<-rbind(ppt.huc12.matrix, ppt.huc12.matrix.sub)
}

nrow(ppt.huc12.matrix)
ppt.huc12.matrix %>% rownames_to_column(var="HUC_ID") %>% 
  mutate(HUC2=substr(HUC_ID,1,4)) %>% 
  group_by(HUC2) %>%
  tally() %>% arrange(desc(n)) 

ppt_hucs<- ppt.huc12.matrix %>% rownames_to_column(var="HUC_ID") %>% pull(HUC_ID)
which(!(hucIDbin$HUC_ID %in% ppt_hucs))

write.csv(ppt.huc12.matrix, '/home/tracidubose/prism_files/ppt_huc12_binary_matrix.csv')

#raster_extract_mean_env should output a matrix 
#that contains 0 where sp is not present and mean environmental values where sp is present
#need to turn that into a long dataframe that excludes 0s

ppt.CS.value<-ppt.huc12.matrix %>% rownames_to_column() %>%
  pivot_longer(-rowname) %>%
  filter(value != 0) %>% 
  rename(scientific_name="name", mean_ppt="value") %>%
  group_by(scientific_name)%>%
  summarize(sd_of_env_mean=sd_2(mean_ppt, na.rm=T))
#check all species have a value (want 104 to print)
nrow(ppt.CS.value %>% filter(!is.na(sd_of_env_mean)))

# problem - if only found in one huc12 SD can't be calculated
# that is why we check above (line 42) for our minimum # of hucs for a sp
write.csv(ppt.CS.value, '/home/tracidubose/rcs_results/HUC_ppt_sd.csv')

test<-ppt.CS.value %>% 
  mutate(ppt_CS_index = (sd_of_env_mean-min(sd_of_env_mean))/
                           (max(sd_of_env_mean)-min(sd_of_env_mean)))


# Extract precipitation data from each species 1km buffer AOO. 
Ppt_Buffer <- data.frame(matrix(NA, nrow = 1, ncol = 144))
colnames(Ppt_Buffer) <- c(substr(FILES_watershed, 1, nchar(FILES_watershed)-10))

for(i in 1:length(FILES_buffer)){
  dat.name <- substr(FILES_buffer[i], 1, nchar(FILES_buffer[i])-4) # Use shapefile name without file extensions.
  spp_BUF <- readOGR(PATH_SHP_1km, dat.name) # Read in dissolved shapefile.
  # Removed buffer b/c of input shapefile.
  BUF_ppt_mean <- raster::extract(mean_ppt, spp_BUF, weights = TRUE, small = TRUE, method = 'simple', df = TRUE) 
  BUF_ppt_mean <- na.omit(BUF_ppt_mean) #There should be no NAs, potential AOO issue with lakes?
  unlist(lapply(BUF_ppt_mean, sd_1))
  Ppt_Buffer[, dat.name] <- sqrt(wtd.var(BUF_ppt_mean$layer, BUF_ppt_mean$weight, na.rm = TRUE, normwt = TRUE))
}

write.csv(Ppt_Buffer, file = "/home/silknets/NFVAP/PRISM/Ppt_Results/Buffer/Ppt_Buffer.csv", row.names=FALSE)
Ppt_SD <- Ppt_Buffer [,-145]
Ppt_SD[2,] <- Watershed_Ppt_SD[1,]
row.names(Ppt_SD) <- c("Buffer", "Watershed")
Ppt_SD <- t(Ppt_SD)
write.csv(Ppt_SD, file = "/home/silknets/NFVAP/PRISM/Ppt_Results/Ppt_SD.csv", row.names=TRUE)


#### MAXIMUM TEMPERATURE CALCS ####

# Set path and load 'raster stack' to variable-specific prism data. 
#options(prism.path = PATH_PRISM_tmax) 
#PRISM_tmax <- prism_stack(ls_prism_data()) # Note: Ensure prism.path is set to correct climate variable

# Mean and Standard Deviation Calculation
#mean_tmax <- raster::calc(PRISM_tmax, fun = mean_1, na.rm = TRUE)
#sd_tmax <- raster::calc(PRISM_tmax, fun = sd_1, na.rm=TRUE)
mean_tmax <- raster("/home/tracidubose/prism_files/Tmax_mean_usa.tif") 
sd_tmax <- raster("/home/tracidubose/prism_files/Tmax_sd_usa.tif")

#Create a dataframe with all 12-digit HUCs and mean precipitation across years. 
# doing this piecemeal using group_hucs; advice: start slow since those are the larger huc groups
tmax.huc12.matrix<-NULL
for(u in group_hucs$HUC2[1:6]){
  hucsub<-hucIDbin %>%
    filter(substr(HUC_ID,1,4)==u)
  
  huc.sf.red<-huc12_red[match(hucsub$HUC_ID, huc12_red$huc12),]
  which(hucsub$HUC_ID != huc.sf.red$huc12)
  
  tmax.huc12.matrix.sub<-extract_mean_env(env.raster=mean_tmax, 
                                         hucs=huc.sf.red, 
                                         hucbinary=hucsub)
  tmax.huc12.matrix<-rbind(tmax.huc12.matrix, tmax.huc12.matrix.sub)
}
#3:55
# good idea to write out ppt.huc12.matrix periodically to save results
write.csv(tmax.huc12.matrix, file = paste0("/home/tracidubose/tmax.incre",
                                           "",#put the increment here so you know
                                           ".csv"), row.names=FALSE)

nrow(tmax.huc12.matrix) #how many hucs have values now?
#which huc groups have you completed?
tmax.huc12.matrix %>% rownames_to_column(var="HUC_ID") %>% 
  mutate(HUC2=substr(HUC_ID,1,4)) %>% 
  group_by(HUC2) %>%
  tally() %>% arrange(desc(n))

#final check
tmax_hucs<- tmax.huc12.matrix %>% rownames_to_column(var="HUC_ID") %>% pull(HUC_ID)
which(!(hucIDbin$HUC_ID %in% tmax_hucs))

write.csv(tmax.huc12.matrix, '/home/tracidubose/prism_files/tmax_huc12_binary_matrix.csv') 
#this is equal to the Watershed_Tmax_SD in Sam's code, 
#file = "/home/silknets/NFVAP/PRISM/Tmax_Results/Watershed/Watershed_Tmax_SD"

#### MINIMUM TEMPERATURE CALCS ####

# Set path and load 'raster stack' to variable-specific prism data. 
#options(prism.path = PATH_PRISM_tmin)
#PRISM_tmin <- prism_stack(ls_prism_data()) # Note: Ensure prism.path is set to correct climate variable

# Mean and Standard Deviation Calculation
#mean_tmin <- raster::calc(PRISM_tmin, fun = mean_1, na.rm = TRUE)
#sd_tmin <- raster::calc(PRISM_tmin, fun = sd_1, na.rm=TRUE)
mean_tmin<-raster("/home/tracidubose/prism_files/Tmin_mean_usa.tif")
sd_tmin <- raster("/home/tracidubose/prism_files/Tmin_sd_usa.tif")

#Create a dataframe with all 12-digit HUCs and mean precipitation across years. 
# doing this piecemeal using group_hucs
tmin.huc12.matrix<-NULL
for(u in group_hucs$HUC2[1:10]){
  hucsub<-hucIDbin %>%
    filter(substr(HUC_ID,1,4)==u)
  
  huc.sf.red<-huc12_red[match(hucsub$HUC_ID, huc12_red$huc12),]
  which(hucsub$HUC_ID != huc.sf.red$huc12)
  
  tmax.huc12.matrix.sub<-extract_mean_env(env.raster=mean_tmin, 
                                          hucs=huc.sf.red, 
                                          hucbinary=hucsub)
  tmin.huc12.matrix<-rbind(tmin.huc12.matrix, tmin.huc12.matrix.sub)
}

# good idea to write out ppt.huc12.matrix periodically to save results
write.csv(tmin.huc12.matrix, file = paste0("/home/tracidubose/tmin.incre",
                                           "",#put the increment here so you know
                                           ".csv"), row.names=FALSE)

nrow(tmin.huc12.matrix) #how many hucs have values now?
#which huc groups have you completed?
tmin.huc12.matrix %>% rownames_to_column(var="HUC_ID") %>% 
  mutate(HUC2=substr(HUC_ID,1,4)) %>% 
  group_by(HUC2) %>%
  tally() %>% arrange(desc(n))

#final check
tmin_hucs<- tmin.huc12.matrix %>% rownames_to_column(var="HUC_ID") %>% pull(HUC_ID)
which(!(hucIDbin$HUC_ID %in% tmin_hucs))

write.csv(tmin.huc12.matrix, '/home/tracidubose/prism_files/tmin_huc12_binary_matrix.csv') 
#this is equal to the Watershed_tmin_SD in Sam's code, 
#file = "/home/silknets/NFVAP/PRISM/Tmin_Results/Watershed/Watershed_tmin_SD"

### stopped here because august max temps, same with jan min temps not stacking?

#### AUGUST MAXIMUM TEMPERATURE CALCS ####

# Set path and load 'raster stack' to variable-specific prism data. 
options(prism.path = PATH_PRISM_tmaxAug) 
PRISM_tmaxAug <- prism_stack(ls_prism_data()) # Note: Ensure prism.path is set to correct climate variable

# Mean and Standard Deviation Calculation
mean_tmaxAug <- calc(PRISM_tmaxAug, fun = mean_1, na.rm = TRUE)
#sd_tmaxAug <- calc(PRISM_tmaxAug, fun = sd_1, na.rm=TRUE)
saveRDS(mean_tmaxAug, file = "/home/silknets/NFVAP/tmaxAug_mean.rds")
#saveRDS(sd_tmaxAug, file = "/home/silknets/NFVAP/tmaxAug_sd.rds")
#mean_tmaxAug <- readRDS(file = "/home/silknets/NFVAP/tmaxAug_mean.rds")
#sd_tmaxAug <- readRDS(file = "/home/silknets/NFVAP/tmaxAug_sd.rds")

# Save Mean and SD Tmax rasters for watersheds. 
#HUC12_tmaxAug_sd <- raster::extract(sd_tmaxAug, HUC12_All, fun= mean ,na.rm = T, weights = TRUE, small = TRUE, method = 'simple', df = TRUE)
#saveRDS(HUC12_tmaxAug_sd, file = "/home/silknets/NFVAP/PRISM/TmaxAug_Results/Watershed/HUC12_tmaxAug_sd")

#write.csv(HUC12_tmaxAug_mean, file = paste0(PATH_NFVAP, "tmaxAug_mean_.csv"), row.names=FALSE)

#Create a dataframe with all 12-digit HUCs and mean precipitation across years. 
HUC12_tmaxAug_mean <- raster::extract(mean_tmaxAug, HUC12_All, fun= mean ,na.rm = T, weights = TRUE, small = TRUE, method = 'simple', df = TRUE)
HUC12_tmaxAug_mean$HUC_ID <- HUC12_All@data$HUC12 
names(HUC12_tmaxAug_mean)[names(HUC12_tmaxAug_mean) == "layer"] <- "mean"
HUC12_tmaxAug_mean <- HUC12_tmaxAug_mean [ , c(3, 2)] 
saveRDS(HUC12_tmaxAug_mean, file = "/home/silknets/NFVAP/PRISM/TmaxAug_Results/Watershed/HUC12_tmaxAug_mean")
HUC12_tmaxAug_mean <- readRDS("/home/silknets/NFVAP/PRISM/TmaxAug_Results/Watershed/HUC12_tmaxAug_mean")

# Make a binary dataframe of occupied watersheds per species.
tmaxAug_Watershed_Occupied <- data.frame(matrix(NA, nrow = 100537, ncol = 145))
colnames(tmaxAug_Watershed_Occupied) <- c("HUC_ID", substr(FILES_watershed, 1, nchar(FILES_watershed)-10))
tmaxAug_Watershed_Occupied$HUC_ID <- HUC12_tmaxAug_mean$HUC_ID

# For loop to populate binary dataframe. 
for(i in 1:length(FILES_watershed)){
  dat.name <- substr(FILES_watershed[i], 1, nchar(FILES_watershed[i])-10)
  Occupied_HUCs <- read.table(file = paste0(PATH_Occupied_HUC12, dat.name, ".txt"), sep = ",", header = T, colClasses = "character")
  tmaxAug_Watershed_Occupied[, dat.name] <- as.integer(tmaxAug_Watershed_Occupied[1:100537,1] %in% Occupied_HUCs$Occupied_HUC)
}

write.csv(tmaxAug_Watershed_Occupied, file = paste0("/home/silknets/NFVAP/data_test_tmaxAug.csv"), row.names=FALSE)

# Use matrix algebra to create a dataframe of SD per species. 
tmaxAug_SD_table <- tmaxAug_Watershed_Occupied [,2:145] * c(HUC12_tmaxAug_mean [,2])

# Build dataframe to be populated with the SD of watershed precipitation per species. 
Watershed_tmaxAug_SD <- data.frame(matrix(NA, nrow = 1, ncol = 144))
colnames(Watershed_tmaxAug_SD) <- c(substr(FILES_watershed, 1, nchar(FILES_watershed)-10))

for(i in 1:length(FILES_watershed)){
  dat.name <- substr(FILES_watershed[i], 1, nchar(FILES_watershed[i])-10)
  Watershed_tmaxAug_SD[, dat.name] <- sd_2(tmaxAug_SD_table[, dat.name], na.rm=TRUE)
}
saveRDS(Watershed_tmaxAug_SD, file = "/home/silknets/NFVAP/PRISM/TmaxAug_Results/Watershed/Watershed_tmaxAug_SD")

# Extract precipitation data from each species 1km buffer AOO. 
tmaxAug_Buffer <- data.frame(matrix(NA, nrow = 1, ncol = 144))
colnames(tmaxAug_Buffer) <- c(substr(FILES_watershed, 1, nchar(FILES_watershed)-10))

for(i in 1:length(FILES_buffer)){
  dat.name <- substr(FILES_buffer[i], 1, nchar(FILES_buffer[i])-4) # Use shapefile name without file extensions.
  spp_BUF <- readOGR(PATH_SHP_1km, dat.name) # Read in dissolved shapefile.
  spp_BUF <- gUnaryUnion(spp_BUF, id=NULL)
  dat.name <- substr(dat.name, 1, nchar(dat.name)-4)
  BUF_tmaxAug_mean <- raster::extract(mean_tmaxAug, spp_BUF, weights = TRUE, small = TRUE, method = 'simple', df = TRUE) # Removed buffer b/c of input shapefile.
  BUF_tmaxAug_mean <- na.omit(BUF_tmaxAug_mean) #There should be no NAs, potential AOO issue with lakes?
  unlist(lapply(BUF_tmaxAug_mean, sd_1))
  tmaxAug_Buffer[, dat.name] <- sqrt(wtd.var(BUF_tmaxAug_mean$layer, BUF_tmaxAug_mean$weight, na.rm = TRUE, normwt = TRUE))
}

write.csv(tmaxAug_Buffer, file = "/home/silknets/NFVAP/PRISM/TmaxAug_Results/Buffer/tmaxAug_Buffer.csv", row.names=FALSE)
tmaxAug_SD <- tmaxAug_Buffer [,-145]
tmaxAug_SD[2,] <- Watershed_tmaxAug_SD[1,]
row.names(tmaxAug_SD) <- c("Buffer", "Watershed")
tmaxAug_SD <- t(tmaxAug_SD)
write.csv(tmaxAug_SD, file = "/home/silknets/NFVAP/PRISM/TmaxAug_Results/tmaxAug_SD.csv", row.names=TRUE)






#### JANUARY MINIMUM TEMPERATURE CALCS ####

# Set path and load 'raster stack' to variable-specific prism data. 
options(prism.path = PATH_PRISM_tminJan) 
PRISM_tminJan <- prism_stack(ls_prism_data()) # Note: Ensure prism.path is set to correct climate variable

# Mean and Standard Deviation Calculation
mean_tminJan <- calc(PRISM_tminJan, fun = mean_1, na.rm = TRUE)
#sd_tminJan <- calc(PRISM_tminJan, fun = sd_1, na.rm=TRUE)
saveRDS(mean_tminJan, file = "/home/silknets/NFVAP/tminJan_mean.rds")
#saveRDS(sd_tminJan, file = "/home/silknets/NFVAP/tminJan_sd.rds")
#mean_tminJan <- readRDS(file = "/home/silknets/NFVAP/tminJan_mean.rds")
#sd_tminJan <- readRDS(file = "/home/silknets/NFVAP/tminJan_sd.rds")

# Save Mean and SD Tmax rasters for watersheds. 
#HUC12_tminJan_sd <- raster::extract(sd_tminJan, HUC12_All, fun= mean ,na.rm = T, weights = TRUE, small = TRUE, method = 'simple', df = TRUE)
#saveRDS(HUC12_tminJan_sd, file = "/home/silknets/NFVAP/PRISM/TminJan_Results/Watershed/HUC12_tminJan_sd")

#write.csv(HUC12_tminJan_mean, file = paste0(PATH_NFVAP, "tminJan_mean_.csv"), row.names=FALSE)

#Create a dataframe with all 12-digit HUCs and mean precipitation across years. 
HUC12_tminJan_mean <- raster::extract(mean_tminJan, HUC12_All, fun= mean ,na.rm = T, weights = TRUE, small = TRUE, method = 'simple', df = TRUE)
HUC12_tminJan_mean$HUC_ID <- HUC12_All@data$HUC12 
names(HUC12_tminJan_mean)[names(HUC12_tminJan_mean) == "layer"] <- "mean"
HUC12_tminJan_mean <- HUC12_tminJan_mean [ , c(3, 2)] 
saveRDS(HUC12_tminJan_mean, file = "/home/silknets/NFVAP/PRISM/TminJan_Results/Watershed/HUC12_tminJan_mean")
HUC12_tminJan_mean <- readRDS("/home/silknets/NFVAP/PRISM/TminJan_Results/Watershed/HUC12_tminJan_mean")

# Make a binary dataframe of occupied watersheds per species.
tminJan_Watershed_Occupied <- data.frame(matrix(NA, nrow = 100537, ncol = 145))
colnames(tminJan_Watershed_Occupied) <- c("HUC_ID", substr(FILES_watershed, 1, nchar(FILES_watershed)-10))
tminJan_Watershed_Occupied$HUC_ID <- HUC12_tminJan_mean$HUC_ID

# For loop to populate binary dataframe. 
for(i in 1:length(FILES_watershed)){
  dat.name <- substr(FILES_watershed[i], 1, nchar(FILES_watershed[i])-10)
  Occupied_HUCs <- read.table(file = paste0(PATH_Occupied_HUC12, dat.name, ".txt"), sep = ",", header = T, colClasses = "character")
  tminJan_Watershed_Occupied[, dat.name] <- as.integer(tminJan_Watershed_Occupied[1:100537,1] %in% Occupied_HUCs$Occupied_HUC)
}

write.csv(tminJan_Watershed_Occupied, file = paste0("/home/silknets/NFVAP/data_test_tminJan.csv"), row.names=FALSE)

# Use matrix algebra to create a dataframe of SD per species. 
tminJan_SD_table <- tminJan_Watershed_Occupied [,2:145] * c(HUC12_tminJan_mean [,2])

# Build dataframe to be populated with the SD of watershed precipitation per species. 
Watershed_tminJan_SD <- data.frame(matrix(NA, nrow = 1, ncol = 144))
colnames(Watershed_tminJan_SD) <- c(substr(FILES_watershed, 1, nchar(FILES_watershed)-10))

for(i in 1:length(FILES_watershed)){
  dat.name <- substr(FILES_watershed[i], 1, nchar(FILES_watershed[i])-10)
  Watershed_tminJan_SD[, dat.name] <- sd_2(tminJan_SD_table[, dat.name], na.rm=TRUE)
}
saveRDS(Watershed_tminJan_SD, file = "/home/silknets/NFVAP/PRISM/TminJan_Results/Watershed/Watershed_tminJan_SD")

# Extract precipitation data from each species 1km buffer AOO. 
tminJan_Buffer <- data.frame(matrix(NA, nrow = 1, ncol = 144))
colnames(tminJan_Buffer) <- c(substr(FILES_watershed, 1, nchar(FILES_watershed)-10))

for(i in 1:length(FILES_buffer)){
  dat.name <- substr(FILES_buffer[i], 1, nchar(FILES_buffer[i])-4) # Use shapefile name without file extensions.
  spp_BUF <- readOGR(PATH_SHP_1km, dat.name) # Read in dissolved shapefile.
  spp_BUF <- gUnaryUnion(spp_BUF, id=NULL)
  dat.name <- substr(dat.name, 1, nchar(dat.name)-4)
  BUF_tminJan_mean <- raster::extract(mean_tminJan, spp_BUF, weights = TRUE, small = TRUE, method = 'simple', df = TRUE) # Removed buffer b/c of input shapefile.
  BUF_tminJan_mean <- na.omit(BUF_tminJan_mean) #There should be no NAs, potential AOO issue with lakes?
  unlist(lapply(BUF_tminJan_mean, sd_1))
  tminJan_Buffer[, dat.name] <- sqrt(wtd.var(BUF_tminJan_mean$layer, BUF_tminJan_mean$weight, na.rm = TRUE, normwt = TRUE))
}

write.csv(tminJan_Buffer, file = "/home/silknets/NFVAP/PRISM/TminJan_Results/Buffer/tminJan_Buffer.csv", row.names=FALSE)
tminJan_SD <- tminJan_Buffer [,-145]
tminJan_SD[2,] <- Watershed_tminJan_SD[1,]
row.names(tminJan_SD) <- c("Buffer", "Watershed")
tminJan_SD <- t(tminJan_SD)
write.csv(tminJan_SD, file = "/home/silknets/NFVAP/PRISM/TminJan_Results/tminJan_SD.csv", row.names=TRUE)






#### CAN ABBY EXPLAIN CLEAN FUNCTION?
# Doing the work of the pivot table 
#clean_dataset <- function(df) {
#  df <- as.data.frame(df)
#  df$year <- substr(df$eventDate, start = 1, stop = 4)
#  df$year <- paste0("_", df$year)
#  df_final <- df %>%
#    gather(key = "prism_name", value = "value", 15:134) %>%
#    rowwise() %>%
#    filter(grepl(pattern = year, x = prism_name))
#  return(df_final)
#}

#### SEEMS LIKE A RULE FOR THE MOVING WINDOW?

# Read in .rds files from output and limit the data to only the instances
# where a point occurs in a particular year, eliminating data that do no follow that rule
#PATH_PRISM_ppt_Results <- "G:/Shared drives/NFVAP - National RCS/R Code/National-RCS/PRISM/ppt/Results"
#filenames <- list.files(path = PATH_PRISM_ppt_Results, full.names = T)
#ldf <- lapply(filenames, readRDS)
#res <- lapply(ldf, clean_dataset)

# Calculate standard deviation on the values for each species and climate variable combination
#### ABBY, ANY SUCCESS WITH THIS CODE?

res <- lapply("G:/Shared drives/NFVAP - National RCS/R Code/National-RCS/PRISM/ppt/Results/Acantharchus pomotis.rds", readRDS)
stdev <- lapply(res[[i]]$value, sd)
stdev <- c()
for (i in 1:length(res)){
  stdev[[i]] <- sd(res[[i]]$value, na.rm = T)# unsure why there are some NAs, not all sp-clim_var have them
  stdev <- list.append(stdev, res[[i]]$name[1])
  stdev <- list.append(stdev, res[[i]]$prism_name[1])
} #need to add in sp and clim_var to list somehow
stdev


###############################

#### CHUNKS BELOW FROM ALB CODE

options(prism.path = "/lustre/projects/css/csas/albenson/prismtmp/ppt")

# Create raster stack for only the years of PRISM data needed
RS <- prism_stack(ls_prism_data()) #raster data
proj4string(RS)<-CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84")

#### ANY REASON WHY BUFFERED .SHP WOULDN"T WORK INSTEAD?
# Run 1km buffers for each point in the species shapefile and extract PRISM data within
# those buffers
bufferSize_ls <- 1000
myextractfun <- function(b){
  result <- raster::extract(RS, finalpoints_sp, buffer = b, fun= mean, na.rm = T, weights = TRUE, small = TRUE, method = 'simple', df = TRUE, sp = TRUE)
}

result <- pbdLapply(bufferSize_ls, myextractfun)
finalize()

# Save result as an RDS
saveRDS(result, file = "ppt_buffers.rds")
