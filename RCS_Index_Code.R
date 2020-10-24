# National Fishes Vulnerability Assessment Project - "RCS_Index_Code.R"
# Revised by Sam Silknetter, 29June2020
# Revised by Traci DuBose, 20Sept2020

# This code calculates Climate Sensitivity (CS) from standard deviation values, 
# merges prior AOO calculations into an RCS output, 
# and calculates Rarity and Climate Sensitivity index (RCS) for each species and range metric. 

# Install and load necessary libraries.
library(tidyverse)

# Import climatic breath for each species and area of occurrence -----
AOOs <- read.csv('rcs_results/AOO HUC12 Output_20201019.csv')
AOO_L48 <- read.csv('rcs_results/AOOs_raw_20201019.csv')
climate_ws_raw <- read.csv('rcs_results/taxa_ws_climate_values20201021.csv')
climate_buf <- read.csv('rcs_results/range.wide.weighted.CS_20201022.csv') 
climate_buf_L48 <- read.csv('rcs_results/buff_climate_sums_20201022.csv')

# Identify taxa we are excluding -----
#   because only 20% of their range (dataframe created within RCS_decision_notes.Rmd)
sp_range_wiL48 <-read.csv('rcs_results/sp_range_within_US.csv') %>% 
  dplyr::select(-X)
NAS_nativetrans<-c('Acris crepitans', 'Anaxyrus boreaus',
                   'Dryophytes cinereus', 'Dryophytes gratiosus',
                   'Dryophytes squirellus', 'Dryophytes wrightorum',
                   'Eleutherodactylus cystignathoides','Lithobates berlandieri',
                   'Lithobates blairi','Lithobates catesbeianus',
                   'Lithobates clamitans','Lithobates grylio',
                   'Lithobates pipiens','Lithobates sphenocephalus',
                   'Pseudacris regilla','Rana aurora',
                   'Rhinella marina')
NAS_exotic<-c('Osteopilus septentrionalis','Xenopus laevis', 'Rhinella marina', 'Lithobates fisheri')
sp_exclude <- sp_range_wiL48 %>% filter(per_inside_us <= .200) %>%
  pull(scientific_name) %>% c(NAS_exotic)
sp_exclude

# Calculate RCS index for entire native range -----
# Scale standard deviations for each climatic variable and grain between 0 and 1. 

# Watershed grain size - climate ---
# how many watersheds for each taxa?
climate_ws_raw %>% filter(value_type=='tmax') %>% 
  group_by(species, value_origin) %>% count()

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

# Buffer grain size - climate ---
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
         SpFac=factor(species, levels = species[order(RCS_WS)])) %>%
  mutate(spatial_extent='entire native NA range')

# Exports RCS data table.
write.csv(RCS_Data, file = paste0("rcs_results/RCS_table_", format(Sys.Date(), "%Y%m%d"),".csv"))

# Calculate RCS index for native range within contiguous US -----
# Scale standard deviations for each climatic variable and grain between 0 and 1. 
# Watershed grain size - climate ---
watershed_CS_L48<-climate_ws_raw %>%
  group_by(species, value_type) %>%
  filter(!(species %in% sp_exclude),
         value_origin == 'huc12') %>%
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

# Buffer grain size - climate ---
buffer_CS_L48 <- climate_buf_L48 %>%
  filter(!(species %in% sp_exclude),
         source=='prism') %>%
  group_by(name) %>% #group_by to scale the climate variables appropriately
  mutate(var_CS_index = (sd_value-min(sd_value))/ #calculating the index
           (max(sd_value)-min(sd_value))) %>% 
  dplyr::select(-mean_value, -sd_value, -sumweight) %>%
  pivot_wider(names_from='name',values_from='var_CS_index') %>% #switching from long to wide dataframe
  group_by(species) %>%
  dplyr::summarise(ppt_buf.CS=mean(ppt, na.rm=T),
                   tmax_buf.CS=mean(tmax, na.rm=T),
                   tmin_buf.CS=mean(tmin, na.rm=T)) %>%
  rowwise() %>%
  mutate(buff_CS=mean(c(ppt_buf.CS, tmax_buf.CS, tmin_buf.CS))) #calculate the mean climate breadth

# Scale AOO values/grain between 0 and 1. 
AOOs_index_L48 <- AOO_L48 %>% 
  filter(!(scientific_name %in% sp_exclude),
         area.type %in% c('usa_l48_alb_1km', 'huc12_WS')) %>%
  dplyr::select(-X, -DFdif) %>%
  filter(!(scientific_name %in% c("Acris blanchardi", 'Lithobates kauffeldi') & area.sqkm ==0)) %>%
  pivot_wider(names_from = area.type, values_from = area.sqkm) %>%
  mutate(WS_AOO_ind = (huc12_WS - min(huc12_WS))/
           (max(huc12_WS) - min(huc12_WS)),
         buff_AOO_ind = (usa_l48_alb_1km - min(usa_l48_alb_1km))/
           (max(usa_l48_alb_1km) - min(usa_l48_alb_1km))) %>%
  rename(species='scientific_name') #for easy joining

# Join the three data frames together (AOOs, WS_CS & buffer_CS)
RCS_Data_L48 <- AOOs_index_L48 %>%
  left_join(buffer_CS_L48, by='species') %>%
  left_join(watershed_CS_L48, by='species') %>%
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
         SpFac=factor(species, levels = species[order(RCS_WS)])) %>%
  mutate(spatial_extent='native L48 range')

# Exports RCS data table.
write.csv(RCS_Data_L48, file = paste0("rcs_results/RCS_table_L48_", format(Sys.Date(), "%Y%m%d"),".csv"))

# Investigate differences between two spatial extents -----
library(ggrepel)
extreme_dif_taxa<-RCS_Data %>%
  bind_rows(RCS_Data_L48) %>% 
  dplyr::select(species, spatial_extent, RCS_WS, RCS_buff) %>%
  pivot_longer(cols=c('RCS_WS', 'RCS_buff')) %>%
  pivot_wider(names_from=spatial_extent, values_from=value) %>%
  group_by(name) %>%
  mutate(rcs_dif=`entire native NA range`-`native L48 range`,
         rank_dif=rank(`entire native NA range`)-rank(`native L48 range`)) %>%
  arrange(rcs_dif) %>% #ungroup()%>%
  filter(species != 'Incilius valliceps') %>%
  dplyr::slice(1:3, 90:93)

RCS_Data %>%
  bind_rows(RCS_Data_L48) %>% 
  dplyr::select(species, spatial_extent, RCS_WS, RCS_buff) %>%
  pivot_longer(cols=c('RCS_WS', 'RCS_buff')) %>%
  pivot_wider(names_from=spatial_extent, values_from=value) %>%
  group_by(name) %>%
  mutate(rcs_dif=`entire native NA range`-`native L48 range`,
         rank_dif=rank(`entire native NA range`)-rank(`native L48 range`)) %>%
  filter(species != 'Incilius valliceps') %>%
  ggplot(aes(x=`entire native NA range`, y=`native L48 range`))+
  geom_polygon(data=data.frame(x=c(0,1,1), y=c(0,0,1)),
               aes(x=x, y=y,fill='less threatened if only consider contiguous US'),  alpha=0.5)+
  geom_point()+
  geom_text_repel(data=extreme_dif_taxa,
                  aes(label=species), size=3,
                  min.segment.length = 0,
                  point.padding=0.6,
                  box.padding = 0.01,
                  fontface='italic')+
  geom_abline(slope=1, intercept=0)+
  facet_wrap(~name, scales='free')+
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(expand=c(0,0))+
  scale_fill_manual('',values='goldenrod3')+
  theme_cowplot()+
  theme(legend.position = 'bottom',
        legend.justification = 'center')
ggsave('rcs_results/figures/sup_diff_btwn_range_rcs.jpg',
       width=6, height=4)

### Plots -----
library(ggplot2); library(cowplot); library(scales)

source('RCS_Conserv_Status.R') #identifies which taxa are on conservation lists
con_plot_df<-RCS_index_con %>%
  mutate(esa=ifelse(`Endangered Species Act`=='Not on ESA', '0', 'ESA'),
         swap=ifelse(`Species of Greatest\nConserv. Need`=='Not in SWAP', '0','SGCN'),
         iucn=ifelse(`IUCN Red List` %in% c('Not listed', 'Least Concern', 'Near Threatened'), '0','IUCN'),
         conlists=paste0(esa, swap, iucn)) %>%
  dplyr::select(species, esa, swap, iucn, conlists) %>%
  mutate(Lissst = case_when(conlists=='ESASGCNIUCN'~'ESA, IUCN, & SGCN',
                                  conlists=="0SGCNIUCN"~'IUCN & SGCN',
                                  conlists=="000"~"None",
                                  T~gsub('0','', conlists)),
         `Listed on:` = factor(Lissst, levels=c("ESA, IUCN, & SGCN", "IUCN & SGCN", 
                                                "SGCN", "IUCN", "None"))) %>%
  filter(!duplicated(.))

unique(con_plot_df$`Listed on:`)
# RCS Dotplot code
RCS_Data %>%
  dplyr::select(species, SpFac, RCS_WS, RCS_buff) %>%
  pivot_longer(cols=c('RCS_WS','RCS_buff')) %>% 
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
  filter(!duplicated(.)) %>%
  #group_by(facet.group)%>% summarize(max(value))
  ggplot()+
  geom_rect(aes(fill=`Listed on:`, xmin=box.min, xmax=box.max, 
                ymin=boy.min, ymax=boy.max))+
  geom_point(aes(x=value, y=SpFac,
                 shape=grain.size), size=2, alpha=0.6)+
  scale_shape_manual('Grain Size', values=c(16,0))+
  scale_fill_manual(values=c(viridis_pal(option='magma', 
                                         begin=.3, end=.75)(3), 
                             'lightgrey'))+
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
        legend.title = element_text(size=12))+
  guides(fill=guide_legend(nrow=2))

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

RCS_Data %>% 
  dplyr::select(species, log10_buff_sqkm, log10_WS_sqkm, 
                WS_AOO_ind, buff_AOO_ind) %>%
  mutate(WS = (log10_WS_sqkm - min(log10_WS_sqkm))/
           (max(log10_WS_sqkm) - min(log10_WS_sqkm)),
         buff = (log10_buff_sqkm - min(log10_buff_sqkm))/
           (max(log10_buff_sqkm) - min(log10_buff_sqkm))) %>%
  dplyr::select(species, WS_AOO_ind, buff_AOO_ind, WS, buff) %>%
  pivot_longer(-species) %>%
  mutate(grain.size=case_when(grepl('WS', name)~'watershed',
                              T~'buffer')) %>%
  ggplot()+
  geom_density(aes(x=value, fill=name), alpha=0.5) +
  facet_wrap(~grain.size)

RCS_Data %>% 
  dplyr::select(species, log10_buff_sqkm, log10_WS_sqkm, 
                WS_AOO_ind, buff_AOO_ind) %>%
  mutate(WS = (log10_WS_sqkm - min(log10_WS_sqkm))/
           (max(log10_WS_sqkm) - min(log10_WS_sqkm)),
         buff = (log10_buff_sqkm - min(log10_buff_sqkm))/
           (max(log10_buff_sqkm) - min(log10_buff_sqkm))) %>%
  dplyr::select(species, WS_AOO_ind, buff_AOO_ind, WS, buff) %>%
  ggplot()+
  geom_point(aes(x=WS_AOO_ind, y=1, color='watershed'), alpha=0.5) +
  geom_point(aes(x=buff_AOO_ind, y=2, color='buff'), alpha=0.5) +
  geom_point(aes(x=WS, y=3, color='watershed'), alpha=0.5) +
  geom_point(aes(x=buff, y=4, color='buff'), alpha=0.5) +
  scale_color_viridis_d()

RCS_Data %>%
  pivot_longer(cols=c("WS_AOO_ind","buff_AOO_ind")) %>%
  ggplot()+geom_density(aes(value, fill=name))

area_order<-RCS_Data %>% arrange(desc(WS_AOO_ind)) %>%
  pull(species)

RCS_Data %>% 
  select(SpFac, species, watershed, buffer) %>%
  mutate(SpFacA=factor(species, levels=area_order),
         filt.rank=as.numeric(SpFacA),
         boy.min=filt.rank-0.4,
         boy.max=filt.rank+0.4) %>%
  left_join(con_plot_df) %>%
  filter(!duplicated(.))%>% 
  ggplot() +
  geom_rect(aes(fill=`Listed on:`, xmin=6, xmax=10, 
                ymin=boy.min, ymax=boy.max))+
  geom_point(aes(x=watershed, y=SpFacA, shape="Watershed"), size=2)+
  geom_point(aes(x=buffer, y=SpFacA, shape="1km Buffer"), size=2)+
  scale_x_log10(name=expression("Area of Occurrence km"^2), expand=c(0,0))+
  scale_y_discrete(name="")+
  scale_shape_manual('Grain Size', values=c(16,0))+
  scale_fill_manual(values=c(viridis_pal(option='magma', 
                                         begin=.3, end=.75)(3), 
                             'lightgrey'))+
  theme_cowplot()+
  theme(panel.grid.major.y = element_line(color="lightgrey"),
        axis.text.y = element_text(size=8, face='italic'),
        axis.text.x = element_text(size=9, angle=20))

ggsave(paste0('rcs_results/figures/area_plot_area_order_jpg_',
              format(Sys.Date(), '%Y%m%d'),'.jpg'), width=6, height=9.5)

# Climate Niche Breadth Plot
climate_order <- RCS_Data %>% arrange(desc(WS_CS)) %>%
  pull(species)
RCS_Data %>% 
  select(SpFac, species, ends_with("CS")) %>%
  pivot_longer(cols=c("ppt_WS.CS", 'tmax_WS.CS','tmin_WS.CS',
                      'ppt_buf.CS','tmax_buf.CS','tmin_buf.CS')) %>%
  mutate(grain.size=case_when(grepl('WS', name)~'watershed',
                              grepl('buf',name)~'buffer'),
         c.var=gsub('_.*','', name),
         SpFacC=factor(species, levels = climate_order)) %>%
  mutate(filt.rank=as.numeric(SpFac),
         boy.min=filt.rank-0.4,
         boy.max=filt.rank+0.4) %>%
  left_join(con_plot_df) %>%
  filter(!duplicated(.))%>% 
  ggplot() +
  geom_rect(aes(fill=`Listed on:`, xmin=-.05, xmax=.00, 
                ymin=boy.min, ymax=boy.max))+
  geom_point(aes(x=value, y=SpFac, color=c.var, shape=c.var), alpha=0.4)+
  geom_point(data=. %>% mutate(grain.size='watershed'),
             aes(x=WS_CS, y=SpFac, color='mean', shape='mean'), size=2)+
  geom_point(data=. %>% mutate(grain.size='buffer'),
             aes(x=buff_CS, y=SpFac, color='mean', shape='mean'),  size=2)+
  scale_x_continuous(name="Climate Breadth Index",
                     expand=c(0,0),
                     labels=c('0','0.25','0.5','0.75','1'))+
  scale_y_discrete(name="")+
  scale_color_manual(name="Climate\nVariable",
                     values = c('black', viridis::viridis_pal()(3)),)+
  scale_shape_manual(name="Climate\nVariable",
                     values=c(8,16,17,15))+
  scale_fill_manual(values=c(viridis_pal(option='magma', 
                                         begin=.3, end=.75)(3), 
                             'lightgrey'))+ 
  theme_cowplot()+
  theme(panel.grid.major.y = element_line(color="lightgrey"),
        axis.text.y = element_text(size=8, face='italic'),
        axis.title.y=element_text(size=-1),
        axis.text.x = element_text(size=9),
        axis.title.x = element_text(size=11),
        strip.background = element_rect(fill=NA),
        strip.text=element_text(size=11))+
  facet_wrap(~grain.size)

ggsave(paste0('rcs_results/figures/CS_plot_rcs_order_jpg_',
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
  facet_wrap(~climate_var, scales="free")+
  scale_fill_manual(values=c('black','white'))+
  theme_cowplot()
