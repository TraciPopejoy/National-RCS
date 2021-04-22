library(tidyverse); library(sf);library(cowplot)

# Set spatial coordinate reference system (CRS).
crs.geo <- st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84")
crs.albers <- st_crs("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs")#2.12kb

huc_bin_long<-read.csv('/home/tracidubose/RCS_Results/HUC12_binary_matrix_20210330.csv') %>%
  pivot_longer(-watershed_ID) %>%
  filter(value!=0) %>%
  mutate(huc12=ifelse(nchar(watershed_ID)==11, paste0(0,watershed_ID),as.character(watershed_ID)),
         taxa=gsub("\\."," ", name))

hucs <- readRDS('/home/tracidubose/huc12_wgs84_sf_02010323.rds') %>% 
  rename(geometry=shape) %>%
  dplyr::select(huc12, geometry, WSAREA, states)
  
# rcs avg value ----
rcs_values<-read.csv("/home/tracidubose/RCS_Results/RCS_table_20210405.csv") %>%
  dplyr::select(-X)
wat_info<-huc_bin_long %>% left_join(rcs_values, by=c('taxa'='species')) %>%
  filter(spatial_extent=='entire native NA range') %>%
  dplyr::select(taxa, huc12, RCS_WS) %>%
  mutate(huc8=substr(huc12, 1,8)) %>%
  group_by(huc8) %>%
  dplyr::summarize(n_taxa=n(),
            avg_rcs=mean(RCS_WS))


huc8<-read_sf(dsn='/home/tracidubose/WBD_National_GDB.gdb',
        layer="WBDHU8")

usal48<-readRDS('/home/tracidubose/usal48_nad83_sf.rds') %>% st_transform(st_crs(huc8))
nat_rcs<-huc8 %>% st_filter(usal48) %>% left_join(wat_info) %>%
  filter(!(states %in% c("MX","CN")),
         huc8 != '04240002') %>%
  ggplot()+
  geom_sf(aes(fill=avg_rcs, color=avg_rcs))+
  theme_bw()+
  scale_color_viridis_c(option='A', name="Average\nRCS value",
                        aesthetics = c('fill','color'),
                        na.value='lightgrey')+
  theme(legend.position = 'bottom')
ggsave('/home/tracidubose/RCS_Results/sup_fig/anuran_rcs_map_0407.jpeg',
       width=3.5, height=3.5)
# species richness ----
huc_count <-huc_bin_long %>% 
  mutate(huc8=substr(huc12, 1,8)) %>% 
  distinct(taxa, huc8, .keep_all = T) %>% 
  count(huc8) 


huc_n<-huc8 %>% left_join(huc_count) %>%
  st_filter(usal48) %>%
  filter(!(states %in% c("MX","CN")),
         huc8 != '04240002')
huc_n %>% filter(states=="CN,MI") %>% ggplot()+geom_sf()+geom_sf_text(aes(label=huc8))

nat_rich<-ggplot(huc_n, mapping=aes(color=n, fill=n))+
  geom_sf()+
  scale_color_viridis_c(option='A', name="Anuran\nRichness",
                        aesthetics = c('fill','color'),
                        na.value='lightgrey')+
  theme_bw()+
  theme(legend.position='bottom')
ggsave('/home/tracidubose/RCS_Results/sup_fig/FigSA_hucdiversity0407.jpeg',
       width=3.5, height=3.5)

fig1<-plot_grid(nat_rich, nat_rcs, labels="AUTO")
ggsave(plot=fig1, '/home/tracidubose/RCS_Results/sup_fig/Fig1_two_maps.jpeg', width=6.5, height=3.5)
