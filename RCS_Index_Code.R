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
NAS_nativetrans<-c('Acris crepitans', 'Anaxyrus boreaus',
                   'Dryophytes cinereus', 'Dryophytes gratiosus',
                   'Dryophytes squirellus', 'Dryophytes wrightorum',
                   'Eleutherodactylus cystignathoides','Lithobates berlandieri',
                   'Lithobates blairi','Lithobates catesbeianus',
                   'Lithobates clamitans','Lithobates grylio',
                   'Lithobates pipiens','Lithobates sphenocephalus',
                   'Pseudacris regilla','Rana aurora',
                   'Rhinella marina',)
NAS_exotic<-c('Osteopilus septentrionalis','Xenopus laevis')
sp_exclude <- sp_range_wiL48 %>% filter(per_inside_us <= .200) %>%
  pull(scientific_name) %>% c(NAS_exotic)
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
         SpFac=factor(species, levels = species[order(RCS_WS)]))

# Exports RCS data table.
write.csv(RCS_Data, file = paste0("rcs_results/RCS_table_", format(Sys.Date(), "%Y%m%d"),".csv"))

### Plots -----
# converting dot plot to ggplot
library(ggplot2); library(cowplot); library(scales)

source('RCS_Conserv_Status.R') #identifies which taxa are on conservation lists
source('investigating statistical options.R') #builds the rcs dataframe
con_plot_df<-rcs %>% dplyr::select(species, RCS_buff,
                                   ESA_bin, State_bin, Int_bin) %>%
  rowwise() %>%
  mutate(connum=sum(ESA_bin, State_bin, Int_bin)) %>%
  left_join(RCS_Data %>% dplyr::select(species, SpFac),
            by='species') %>%
  arrange(RCS_buff) %>% ungroup() %>%
  filter(!duplicated(.))

# RCS Dotplot code
RCS_Data %>%
  dplyr::select(SpFac, RCS_WS, RCS_buff) %>%
  pivot_longer(-SpFac) %>% 
  mutate(grain.size=recode(name, RCS_buff='1km Buffers',
                           RCS_WS='Watersheds'),
         filt.rank=as.numeric(SpFac),
         box.min=ifelse(filt.rank > 49, 1.01, .85),
         box.max=ifelse(filt.rank > 49, 1.03, .899),
         boy.min=ifelse(filt.rank > 49, filt.rank-0.4-49, filt.rank-0.4),
         boy.max=ifelse(filt.rank > 49, filt.rank+0.4-49, filt.rank+0.4),
         facet.group=factor(case_when(filt.rank > 49~'vulnerable',
                               T~'not as vulnerable'),
                               levels=c('vulnerable','not as vulnerable'))) %>%
  left_join(con_plot_df) %>%
  #group_by(facet.group)%>% summarize(max(value))
  ggplot()+
  geom_rect(aes(fill=as.factor(connum), xmin=box.min, xmax=box.max, 
                ymin=boy.min, ymax=boy.max))+
  geom_point(aes(x=value, y=SpFac,
                 shape=grain.size), size=2, alpha=0.6)+
  scale_shape_manual('Grain Size', values=c(16,0))+
  scale_fill_manual('# Conservation Lists',
                    values=c("black","#440154FF","#2A788EFF","#7AD151FF"))+
  scale_x_reverse(name="RCS Index", expand=c(0,0), 
                  breaks=c(0, .25, .5, .7, .8, .9, 1))+
  scale_y_discrete(name="")+
  theme_cowplot()  +
  facet_wrap(~facet.group, scales = 'free') +
  theme(panel.grid.major.y = element_line(color="lightgrey"),
        axis.text.y = element_text(size=8, face='italic'),
        axis.title.y=element_text(size=-1),
        axis.text.x = element_text(size=10, angle=30, hjust=.9),
        axis.title.x = element_text(size=10),
        strip.background = element_rect(fill=NA),
        strip.text=element_text(size=11),
        legend.position='bottom',
        legend.box="vertical",
        legend.title = element_text(size=12))

ggsave(paste0('rcs_results/figures/RCS_jpg_new_', 
              format(Sys.Date(), '%Y%m%d'),'.jpg'), width=6, height=7)

count_rcs<-RCS_Data %>%
  left_join(read.csv('data/anuran_occ_all20200708.csv', stringsAsFactors = F) %>%
              count(final.taxa),
            by=c('species'='final.taxa')) %>%
  dplyr::select(SpFac, RCS_WS, RCS_buff, n) %>%
  pivot_longer(cols = c('RCS_WS', 'RCS_buff'))
 
ggplot(data=count_rcs)+
  geom_point(aes(x=value, y=n))+
  scale_y_log10()+
  geom_smooth(aes(x=value, y=n), method='lm')+
  facet_wrap(~name)

summary(lm(value~n+name, data=count_rcs))

ggplot()+
  geom_point(data=RCS_Data, 
             aes(x=log10_WS_sqkm, y=SpFac, shape="Watershed"), size=2)+
  geom_point(data=RCS_Data, 
             aes(x=log10_buff_sqkm, y=SpFac, shape="1km Buffer"), size=2)+
  scale_x_continuous(name="log(Area of Occurrence)")+
  scale_y_discrete(name="")+
  scale_shape_manual('Grain Size', values=c(16,0))+
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
