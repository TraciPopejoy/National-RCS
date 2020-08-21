# converting dot plot to ggplot

library(tidyverse); library(ggplot2); library(cowplot);library(scales)

# Read in RCS data
RCS_Data <- read.csv("rcs_results/RCS_table_20200726.csv") %>%
  dplyr::select(-X)
  #adding a species name factor that is ordered based on RCS rank

# important variables I might want to plot: 
# AOO_WS_adj, CS_WS_adj, AOO_BUF_adj, CS_BUF_adj, RCS_WS, RCS_BUF

#adding a species name factor that is ordered based on RCS rank
species_ordered_RCS <- RCS_Data[order(RCS_Data$RCS_WS, decreasing = F),]$scientific_name
species_ordered_RCS <- gsub("\\.", " ", species_ordered_RCS) #replacing period with space
RCS_Data<-RCS_Data %>% mutate(SpFac=factor(gsub("\\."," ", scientific_name),
                                     levels=species_ordered_RCS)) %>%
  dplyr::select(scientific_name, SpFac, everything()) #reordering for ease of checking df

# RCS Dotplot code
ggplot()+
  geom_point(data=RCS_Data, aes(x=RCS_WS, y=SpFac))+
  scale_x_reverse(name="RCS Index")+
  scale_y_discrete(name="")+
  theme_cowplot()  +
  theme(panel.grid.major.y = element_line(color="lightgrey"),
        axis.text.y = element_text(size=10),
        axis.title.x = element_text(size=15))
ggsave('rcs_results/figures/RCS_jpg.jpg', width=7, height=15)

ggplot()+
  geom_point(data=RCS_Data, 
             aes(x=logHUC12, y=SpFac,color="HUC 12"), size=2)+
  geom_point(data=RCS_Data, 
             aes(x=logBuffer1km, y=SpFac, color="1km buffer"), size=2)+
  scale_x_continuous(name="log(Area of Occurrence)")+
  scale_y_discrete(name="")+
  scale_color_viridis_d(name="Grain Size",option="plasma",
                        end=.6)+
  theme_cowplot()+
  theme(panel.grid.major.y = element_line(color="lightgrey"),
        axis.text.y = element_text(size=10),
        axis.title.x = element_text(size=15),
        legend.position = 'top')
ggsave('rcs_results/figures/area_plot_jpg.jpg', width=7, height=15)

RCS_climate_plot<-RCS_Data %>% 
  select(scientific_name, SpFac,ends_with("CS")) %>%
  pivot_longer(cols=c("ppt_CS", 'tmax_CS','tmin_CS','tmean_CS'))
ggplot()+
  geom_point(data=RCS_climate_plot, 
             aes(x=value, y=SpFac,color=name, shape=name), alpha=0.7)+
  geom_point(data=RCS_Data, 
             aes(x=WS_CS, y=SpFac, color='mean',shape='mean'), size=2)+
  scale_x_continuous(name="Climate Breadth Index")+
  scale_y_discrete(name="")+
  scale_color_manual(name="Climate\nVariable",
                     values = c('black', viridis_pal()(4)),
                     labels=c('mean CS','ppt CS','tmax CS',
                              'tmean CS', 'tmin CS'))+
  scale_shape_manual(name="Climate\nVariable",
                     values=c(1,16,16,16,16),
                     labels=c('mean CS','ppt CS','tmax CS',
                              'tmean CS', 'tmin CS'))+
  theme_cowplot()+
  theme(panel.grid.major.y = element_line(color="lightgrey"),
        axis.text.y = element_text(size=10),
        axis.title.x = element_text(size=15))
ggsave('rcs_results/figures/CS_plot_jpg.jpg', width=7, height=15)
