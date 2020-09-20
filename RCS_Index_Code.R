# National Fishes Vulnerability Assessment Project - "RCS_Index_Code.R"
# Revised by Sam Silknetter, 29June2020
# Revised by Traci DuBose, 20Sept2020

# This code calculates Climate Sensitivity (CS) from standard deviation values, 
# merges prior AOO calculations into an RCS output, 
# and calculates Rarity and Climate Sensitivity index (RCS) for each species and range metric. 

# Install and load necessary libraries.
library(tidyverse)

# Import climatic breath for each species and area of occurrence
AOOs <- read.csv('rcs_results/AOO HUC12 Output_20200920.csv')
climate_ws_raw <- read.csv('rcs_results/taxa_ws_climate_values20200916.csv')
climate_buf <- read.csv('rcs_results/range.wide.weighted.CS_20200920.csv') 

# Scale standard deviations for each climatic variable and grain (n=2) between 0 and 1. ###
# Watershed grain size - climate
watershed_CS<-climate_ws_raw %>%
  group_by(species, value_type) %>%
  dplyr::summarize(mean_climate=mean(value), sd_of_env_mean=sd(value)) %>%
  group_by(value_type) %>% #group_by to scale the climate variables appropriately
  mutate(var_CS_index = (sd_of_env_mean-min(sd_of_env_mean))/ #calculating the index
           (max(sd_of_env_mean)-min(sd_of_env_mean))) %>% 
  pivot_wider(names_from='value_type',values_from='var_CS_index') %>% #switching from long to wide dataframe
  dplyr::select(-mean_climate, -sd_of_env_mean) %>%
  group_by(species) %>%
  dplyr::summarise(ppt_WS.CS=mean(`annual ppt`, na.rm=T),
                   tmax_WS.CS=mean(tmax, na.rm=T),
                   tmin_WS.CS=mean(tmin, na.rm=T)) %>%
  rowwise() %>%
  mutate(WS_CS=mean(c(ppt_WS.CS, tmax_WS.CS, tmin_WS.CS))) #calculate the mean climate breadth
# Buffer grain size - climate
buffer_CS <- climate_buf %>% 
  group_by(name) %>% #group_by to scale the climate variables appropriately
  mutate(var_CS_index = (sd_value-min(sd_value))/ #calculating the index
           (max(sd_value)-min(sd_value))) %>% 
  pivot_wider(names_from='name',values_from='var_CS_index')%>% #switching from long to wide dataframe
  dplyr::select(-mean_value, -sd_value) %>%
  group_by(species) %>%
  dplyr::summarise(ppt_buf.CS=mean(ppt, na.rm=T),
                   tmax_buf.CS=mean(tmax, na.rm=T),
                   tmin_buf.CS=mean(tmin, na.rm=T)) %>%
  rowwise() %>%
  mutate(buff_CS=mean(c(ppt_buf.CS, tmax_buf.CS, tmin_buf.CS))) #calculate the mean climate breadth

# Scale AOO values/grain between 0 and 1. 
AOOs_index<-AOOs %>% 
  dplyr::select(-X, -mean_rank, -sd_rank) %>%
  mutate(WS_AOO_ind = (watershed - min(watershed))/
           (max(watershed) - min(watershed)),
         buff_AOO_ind = (buffer - min(buffer))/
           (max(buffer) - min(buffer))) %>%
  rename(species='scientific_name') #for easy joining

# Join the three data frames together (AOOs, WS_CS & buffer_CS)
RCS_Data <- AOOs_index %>%
  left_join(buffer_CS, by='species') %>%
  left_join(watershed_CS, by='species') %>%
  # Subtract the scaled AOO values from 1 
  # (so that low values = commonness, high values = vulnerability.
  mutate(AOO_WS_adj = 1 - WS_AOO_ind,
         CS_WS_adj = 1 - WS_CS,
         AOO_BUF_adj = 1 - buff_AOO_ind,
         CS_BUF_adj = 1 - buff_CS,
  # Calculate Relative Climate Sensitivity Index by taking average of "1-AOO" and "1-CS".
  #  Where large values of RCS indicate species with small AOO and low range of climate variables.
         RCS_WS = (AOO_WS_adj + CS_WS_adj)/2,
         RCS_buff = (AOO_BUF_adj + CS_BUF_adj)/2)

# Exports RCS data table.
write.csv(RCS_Data, file = paste0("rcs_results/RCS_table_", format(Sys.Date(), "%Y%m%d"),".csv"))


### Plots -----
# converting dot plot to ggplot
library(tidyverse); library(ggplot2); library(cowplot);library(scales)

# Read in RCS data
RCS_Data <- read.csv("rcs_results/RCS_table_20200920.csv") %>%
  dplyr::select(-X)

# important variables I might want to plot: 
# AOO_WS_adj, CS_WS_adj, AOO_BUF_adj, CS_BUF_adj, RCS_WS, RCS_BUF

#adding a species name factor that is ordered based on RCS rank
species_ordered_RCS <- RCS_Data[order(RCS_Data$RCS_WS, decreasing = F),]$species
RCS_Data<-RCS_Data %>% mutate(SpFac=factor(species,
                                           levels=species_ordered_RCS)) %>%
  dplyr::select(species, SpFac, everything()) #reordering for ease of checking df

# RCS Dotplot code
ggplot()+
  geom_point(data=RCS_Data, aes(x=RCS_WS, y=SpFac))+
  scale_x_reverse(name="RCS Index")+
  scale_y_discrete(name="")+
  theme_cowplot()  +
  theme(panel.grid.major.y = element_line(color="lightgrey"),
        axis.text.y = element_text(size=10),
        axis.title.x = element_text(size=15))
ggsave(paste0('rcs_results/figures/RCS_jpg_', 
              format(Sys.Date(), '%Y%m%d'),'.jpg'), width=7, height=15)

ggplot()+
  geom_point(data=RCS_Data, 
             aes(x=log10_WS_sqkm, y=SpFac, color="Watersheds"), size=2)+
  geom_point(data=RCS_Data, 
             aes(x=log10_buff_sqkm, y=SpFac, color="1km buffer"), size=2)+
  scale_x_continuous(name="log(Area of Occurrence)")+
  scale_y_discrete(name="")+
  scale_color_viridis_d(name="Grain Size",option="plasma",
                        end=.6)+
  theme_cowplot()+
  theme(panel.grid.major.y = element_line(color="lightgrey"),
        axis.text.y = element_text(size=10),
        axis.title.x = element_text(size=15),
        legend.position = 'top')
ggsave(paste0('rcs_results/figures/area_plot_jpg_',
              format(Sys.Date(), '%Y%m%d'),'.jpg'), width=7, height=15)

RCS_climate_plot<-RCS_Data %>% 
  select(SpFac, ends_with("CS")) %>%
  pivot_longer(cols=c("ppt_WS.CS", 'tmax_WS.CS','tmin_WS.CS',
                      'ppt_buf.CS','tmax_buf.CS','tmin_buf.CS')) %>%
  mutate(grain.size=case_when(grepl('WS', name)~'watershed',
                              grepl('buf',name)~'buffer'))
ggplot()+
  geom_point(data=RCS_climate_plot, 
             aes(x=value, y=SpFac,color=name, shape=name), alpha=0.7)+
  geom_point(data=RCS_Data, 
             aes(x=WS_CS, y=SpFac, color='mean',), size=2)+
  geom_point(data=RCS_Data, 
             aes(x=buff_CS, y=SpFac, color='mean'), size=2)+
  scale_x_continuous(name="Climate Breadth Index")+
  scale_y_discrete(name="")+
  #scale_color_manual(name="Climate\nVariable",
  #                   values = c('black', viridis_pal()(4)),
  #                   labels=c('mean CS','ppt CS','tmax CS',
  #                            'tmean CS', 'tmin CS'))+
  #scale_shape_manual(name="Climate\nVariable",
  #                   values=c(1,16,16,16,16),
  #                   labels=c('mean CS','ppt CS','tmax CS',
  #                            'tmean CS', 'tmin CS'))+
  theme_cowplot()+
  #theme(panel.grid.major.y = element_line(color="lightgrey"),
  #      axis.text.y = element_text(size=10),
  #      axis.title.x = element_text(size=15))
  facet_wrap(~grain.size)
ggsave(paste0('rcs_results/figures/CS_plot_jpg_',
              format(Sys.Date(), '%Y%m%d'),'.jpg'), width=7, height=15)
