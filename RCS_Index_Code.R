# National Fishes Vulnerability Assessment Project - "RCS_Index_Code.R"
# Revised by Sam Silknetter, 29June2020
# Revised by Traci DuBose, 24July2020

# This code calculates Climate Sensitivity (CS) from standard deviation values, merges prior AOO calculations into
#  an RCS output, and calculates Rarity and Climate Sensitivity index (RCS) for each species and range metric. 

# Install and load necessary libraries.
library(tidyverse)

# Import standard deviations for climatic data (indicate climate breath for each species).
ppt_SD <- read.csv(file = "data/HUC_ppt_sd.csv") %>% select(-X)
Tmax_SD <- read.csv(file = "data/HUC_tmax_sd.csv") %>% select(-X)
Tmin_SD <- read.csv(file = "data/HUC_tmin_sd.csv") %>% select(-X)
Tmean_SD <- read.csv(file = "data/HUC_tmean_sd.csv") %>% select(-X) 

# Scale standard deviations for each climatic variable (n=5) and grain (n=2) between 0 and 1. ###
#HUC12 grain size
huc_cs<-bind_rows(ppt_SD %>% mutate(climate_variable="ppt_CS"), #joining each climate variable sd
          Tmax_SD %>% mutate(climate_variable="tmax_CS"),
          Tmin_SD %>% mutate(climate_variable="tmin_CS"),
          Tmean_SD %>% mutate(climate_variable="tmean_CS")) %>%
  group_by(climate_variable) %>% #group_by to scale the climate variables appropriately
  mutate(var_CS_index = (sd_of_env_mean-min(sd_of_env_mean))/ #calculating the index
                        (max(sd_of_env_mean)-min(sd_of_env_mean))) %>%
  dplyr::select(-sd_of_env_mean) %>% # removing the climate sd raw values
  pivot_wider(names_from='climate_variable',values_from='var_CS_index') %>% #switching from long to wide dataframe
  rowwise() %>%
  mutate(WS_CS=mean(c(ppt_CS,tmax_CS, tmin_CS, tmean_CS))) #calculate the mean climate breadth

# Read in AOO dataframe as RCS Data.
RCS_Data <- read.csv(file = "rcs_results/AOO HUC12 Output_20200721.csv") %>%
  dplyr::select(-X, -sd_Rank, -DFdif, -cv, -mean_Rank, -range, -StandardizedMean) %>%
  mutate(scientific_name=sub(" ",".",scientific_name)) %>%
  full_join(huc_cs, by='scientific_name') #%>%
  #full_join(buff_cs, by='scientific_name')

# Scale AOO values/grain between 0 and 1.
RCS_Data$WS_AOO_scaled <- ((RCS_Data$huc12_area - min(RCS_Data$huc12_area))
  /(max(RCS_Data$huc12_area) - min(RCS_Data$huc12_area)))
RCS_Data$BUF_AOO_scaled <- ((RCS_Data$dissolved_buffer_1km - min(RCS_Data$dissolved_buffer_1km))
  /(max(RCS_Data$dissolved_buffer_1km) - min(RCS_Data$dissolved_buffer_1km)))

# Subtract the scaled AOO values from 1 (so that low values = commonness, high values = vulnerability.
RCS_Data$AOO_WS_adj <- (1 - RCS_Data$WS_AOO_scaled)
RCS_Data$CS_WS_adj <- (1 - RCS_Data$WS_CS)
RCS_Data$AOO_BUF_adj <- (1 - RCS_Data$BUF_AOO_scaled)
RCS_Data$CS_BUF_adj <- (1 - RCS_Data$BUF_CS)

# Calculate Relative Climate Sensitivity Index by taking average of "1-AOO" and "1-CS".
#  Where large values of RCS indicate species with small AOO and low range of climate variables.
RCS_Data$RCS_WS <- ((RCS_Data$AOO_WS_adj + RCS_Data$CS_WS_adj)/2)
RCS_Data$RCS_BUF <- ((RCS_Data$AOO_BUF_adj + RCS_Data$CS_BUF_adj)/2)

# Exports RCS data table.
write.csv(RCS_Data, file = paste0("rcs_results/RCS_table_", format(Sys.Date(), "%Y%m%d"),".csv"))
