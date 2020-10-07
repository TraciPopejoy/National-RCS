# National Fishes Vulnerability Assessment Project - "RCS_Index_Code.R"
# Revised by Sam Silknetter, 29June2020
# Revised by Traci DuBose, 20Sept2020

# This code calculates Climate Sensitivity (CS) from standard deviation values, 
# merges prior AOO calculations into an RCS output, 
# and calculates Rarity and Climate Sensitivity index (RCS) for each species and range metric. 

# Install and load necessary libraries.
library(tidyverse)

# Import climatic breath for each species and area of occurrence
AOOs <- read.csv('rcs_results/AOO HUC12 Output_20200929.csv')
climate_ws_raw <- read.csv('rcs_results/taxa_ws_climate_values20200929.csv')
climate_buf <- read.csv('rcs_results/range.wide.weighted.CS_20200929.csv') 

#identify taxa we are excluding because only 20% of their range?
# dataframe created from code within RCS_decision_notes.Rmd
sp_range_wiL48 <-read.csv('rcs_results/sp_range_within_US.csv') %>% dplyr::select(-X)
sp_exclude <- sp_range_wiL48 %>% filter(per_inside_us <= .200) %>%
  pull(scientific_name)

# Scale standard deviations for each climatic variable and grain (n=2) between 0 and 1. ###
# Watershed grain size - climate
watershed_CS<-climate_ws_raw %>%
  group_by(species, value_type) %>%
  filter(!(species %in% sp_exclude)) %>%
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
  filter(!(species %in% sp_exclude)) %>%
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
  filter(!(scientific_name %in% sp_exclude)) %>%
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
         RCS_buff = (AOO_BUF_adj + CS_BUF_adj)/2,
  #create a factor for ordering plots
         SpFac=factor(species, levels = species[order(RCS_buff)]))

# Exports RCS data table.
write.csv(RCS_Data, file = paste0("rcs_results/RCS_table_", format(Sys.Date(), "%Y%m%d"),".csv"))

### Plots -----
# converting dot plot to ggplot
library(ggplot2); library(cowplot); library(scales)

# RCS Dotplot code
RCS_Data %>%
  dplyr::select(SpFac, RCS_WS, RCS_buff) %>%
  pivot_longer(-SpFac) %>% 
  mutate(grain.size=recode(name, RCS_buff='1km Buffers',
                           RCS_WS='Watersheds'),
         filt.rank=as.numeric(SpFac))%>%
  #filter(filt.rank > 52) %>%
  ggplot()+
  geom_point(aes(x=value, y=SpFac))+
  scale_x_reverse(name="RCS Index")+
  scale_y_discrete(name="")+
  theme_cowplot()  +
  facet_wrap(~grain.size) +
  theme(panel.grid.major.y = element_line(color="lightgrey"),
        axis.text.y = element_text(size=8, face='italic'),
        axis.title.y=element_text(size=-1),
        axis.text.x = element_text(size=10),
        axis.title.x = element_text(size=11),
        strip.background = element_rect(fill=NA),
        strip.text=element_text(size=11))

ggsave(paste0('rcs_results/figures/RCS_jpg_new_', 
              format(Sys.Date(), '%Y%m%d'),'.jpg'), width=6, height=9)

ggplot()+
  geom_point(data=RCS_Data, 
             aes(x=log10_WS_sqkm, y=SpFac, shape="Watershed"), size=2)+
  geom_point(data=RCS_Data, 
             aes(x=log10_buff_sqkm, y=SpFac, shape="1km Buffer"), size=2)+
  scale_x_continuous(name="log(Area of Occurrence)")+
  scale_y_discrete(name="")+
  scale_shape_discrete(name='Grain Size')+
  theme_cowplot()+
  theme(panel.grid.major.y = element_line(color="lightgrey"),
        axis.text.y = element_text(size=10, face='italic'),
        axis.title.x = element_text(size=15),
        legend.position = 'top')
ggsave(paste0('rcs_results/figures/area_plot_jpg_',
              format(Sys.Date(), '%Y%m%d'),'.jpg'), width=6, height=9.5)

RCS_Data %>% 
  select(SpFac, ends_with("CS")) %>%
  pivot_longer(cols=c("ppt_WS.CS", 'tmax_WS.CS','tmin_WS.CS',
                      'ppt_buf.CS','tmax_buf.CS','tmin_buf.CS')) %>%
  mutate(grain.size=case_when(grepl('WS', name)~'watershed',
                              grepl('buf',name)~'buffer'),
         c.var=gsub('_.*','', name)) %>%
  ggplot()+
  geom_point(aes(x=value, y=SpFac,color=c.var, shape=c.var), alpha=0.7)+
  geom_point(data=RCS_Data %>% mutate(grain.size='watershed'),
             aes(x=WS_CS, y=SpFac, color='mean', shape='mean'), size=2)+
  geom_point(data=RCS_Data %>% mutate(grain.size='buffer'),
             aes(x=buff_CS, y=SpFac, color='mean', shape='mean'),  size=2)+
  scale_x_continuous(name="Climate Breadth Index")+
  scale_y_discrete(name="")+
  scale_color_manual(name="Climate\nVariable",
                     values = c('black', viridis_pal()(3)),)+
  scale_shape_manual(name="Climate\nVariable",
                     values=c(8,16,17,15))+
  theme_cowplot()+
  theme(panel.grid.major.y = element_line(color="lightgrey"),
        axis.text.y = element_text(size=8, face='italic'),
        axis.title.y=element_text(size=-1),
        axis.text.x = element_text(size=10),
        axis.title.x = element_text(size=11),
        strip.background = element_rect(fill=NA),
        strip.text=element_text(size=11))+
  facet_wrap(~grain.size)
ggsave(paste0('rcs_results/figures/CS_plot_jpg_',
              format(Sys.Date(), '%Y%m%d'),'.jpg'), width=7, height=9.5)


### investigating difference between buffer and watershed rank rcs ----
names(RCS_Data)
big_sp_shifts_area<- RCS_Data %>%
  summarize(species=species,
            ws_cs_r = rank(WS_CS),
         bf_cs_r = rank(buff_CS),
         aoo_rank_dif = rank_WS-rank_buff,
         cs_rank_dif = ws_cs_r-bf_cs_r) %>%
  as_tibble() %>%
  arrange(aoo_rank_dif) %>%
  #left_join(sp_range_wiL48 %>%
  #            dplyr::select(scientific_name, per_inside_us), 
  #          by=c('species'='scientific_name')) %>%
  filter(aoo_rank_dif > 9 | aoo_rank_dif < -9 ) %>% pull(species)
big_sp_shifts_cs<-RCS_Data %>%
  summarize(species=species,
            ws_cs_r = rank(WS_CS),
            bf_cs_r = rank(buff_CS),
            aoo_rank_dif = rank_WS-rank_buff,
            cs_rank_dif = ws_cs_r-bf_cs_r) %>%
  as_tibble() %>%
  arrange(aoo_rank_dif) %>%
  #left_join(sp_range_wiL48 %>%
  #            dplyr::select(scientific_name, per_inside_us), 
  #          by=c('species'='scientific_name')) %>%
  #ggplot()+geom_histogram(aes(x=cs_rank_dif))
  filter(cs_rank_dif > 20 | cs_rank_dif < -20 ) %>% pull(species)
RCS_Data %>%
  summarize(species=species,
         r_rcs_ws=rank(RCS_WS),
         r_rcs_bf=rank(RCS_buff),
         RCS_dif=r_rcs_ws-r_rcs_bf,
            ws_cs_r = rank(WS_CS),
            bf_cs_r = rank(buff_CS),
            grain_dif = ws_cs_r - bf_cs_r,
            aoo_rank_dif = rank_WS-rank_buff,
            cs_rank_dif = ws_cs_r-bf_cs_r) %>%
  as_tibble() %>% dplyr::select(species, ws_cs_r, bf_cs_r, cs_rank_dif)
#  filter(species %in% unique(c(big_sp_shifts_cs, big_sp_shifts_area))) %>%

 %>%
  summarize(mean(abs(r_rcs_ws-r_rcs_bf)),
            range(r_rcs_ws-r_rcs_bf),
            mean(abs(aoo_rank_dif)),
            range(aoo_rank_dif),
            mean(abs(cs_rank_dif)),
            range(cs_rank_dif))


RCS_Data %>% dplyr::select(species, RCS_WS, RCS_buff) %>% 
  mutate(rank.ws=rank(RCS_WS),
         rank.buff=rank(RCS_buff),
         difff=rank.ws-rank.buff) %>%
  slice(1:7)

climate_raw_values<-bind_rows(climate_buf %>% mutate(grain.size="BUFF"),
                              climate_ws_raw %>% group_by(species, value_type) %>%
                                dplyr::summarize(mean_climate=mean(value), 
                                                 sd_of_env_mean=sd(value)) %>%
                                rename(name='value_type',
                                       mean_value='mean_climate',
                                       sd_value='sd_of_env_mean') %>%
                                mutate(grain.size='WS')) %>%
  mutate(climate_var=recode(name, 'ppt'='annual ppt'))

ggplot(climate_raw_values)+
  geom_density(aes(x=sd_value, fill=grain.size), alpha=.4)+
  facet_wrap(~climate_var, scales="free")
