library(tidyverse); library(sf); library(Hmisc); library(cowplot)

# Set spatial coordinate reference system (CRS).
crs.geo <- st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84")
crs.albers <- st_crs("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs")#2.12kb
PATH_iucn_ranges <-'/home/tracidubose/ANURA/'
iucn_anuran<-read_sf(PATH_iucn_ranges, 'ANURA')  %>% #IUCN maps for native ranges
  st_transform(crs.albers)

# trying to recreate USGS regions -----
Northwest <- c('washington', 'oregon', 'idaho')
Southeast <- c('iowa', 'missouri', 'arkansas', 'oklahoma', 'texas',
               'louisiana', 'mississippi', 'alabama', 'georgia', 'florida',
               'south carolina', 'north carolina', 'tennessee')
Midcontintent <- c('ohio', 'indiana', 'illinois', 'michigan', 
                   'wisconsin','minnesota', 'kansas',
                   'south dakota', 'north dakota','nebraska', 'montana')
RockyMountain <- c('wyoming', 'colorado', 'utah', 'new mexico')
Southwest <- c('arizona', 'california', 'nevada')
Northeast <- c('maine', 'new york', 'vermont', 'new hampshire', 'massachusetts',
               'delaware', 'connecticut', 'district of columbia', 'kentucky',
               'maryland', 'new jersey', 'pennsylvania', 'rhode island', 
               'west virginia','virginia')


library(maps)
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE)) %>% 
  mutate(usgs_region=case_when(ID %in% Northwest ~ 'Northwest',
                               ID %in% Northeast ~ 'Northeast',
                               ID %in% Southwest ~ 'Southwest',
                               ID %in% Southeast ~ 'Southeast',
                               ID %in% Midcontintent ~ 'Midcontinent',
                               ID %in% RockyMountain ~ 'Rocky Mountain')) %>%
  st_transform(crs.albers)
ggplot()+geom_sf(data=states, aes(fill=usgs_region))

#split regions apart
for(i in unique(states$usgs_region)){
  assign(paste0(i, '_sf'),
         states %>% filter(usgs_region == i) %>%
           st_buffer(0) %>%
           st_union())
}
usgs_reg_sf<-paste0(unique(states$usgs_region), '_sf')
ggplot()+geom_sf(data=Southeast_sf)

# for-loop to create 4 different shapefiles for each anuran taxa -----
AOOs<-NULL
short_run<-c("Lithobates yavapaiensis", "Lithobates onca", "Dryophytes wrightorum","Dryophytes versicolor")
PATH_regional_spp_shp = '/home/tracidubose/rcs_results/regional_species_shapefiles/'
pnw_taxa<-c("Anaxyrus woodhousii", "Ascaphus montanus", "Rana luteiventris", 
            "Lithobates pipiens", "Spea intermontana", "Rana pretiosa",
            "Rana boylii", "Rana cascadae", "Rana aurora", "Ascaphus truei",
            "Anaxyrus boreas", "Pseudacris regilla")
entire_range_rcs<-read.csv('/home/tracidubose/RCS_table_20201028.csv')
anuran.taxa <- entire_range_rcs %>% pull(species)
PATH_FocalSpecies_OccurrenceData <- "/home/tracidubose/RCS_Anuran_input/Anuran Occurrence Data/"

# Calculate Area of Occupancy ----
reg_aoo<-NULL
for(y in usgs_reg_sf){
  for(i in anuran.taxa){
    reg_out<-find_area(taxa=i, spatial.file1=y,
                       output.path=paste0(PATH_regional_spp_shp, sub('_sf', '', y), '/'),
                       watershed.name.column = NA)
    reg_aoo<-bind_rows(reg_aoo,reg_out)
    print(paste(i, y))
    rm(reg_out)
  }
}
reg_aoo %>% 
  pivot_wider(names_from=area.type, values_from = area.sqkm) %>% 
  mutate(rowS=rowSums(.[-1])) %>%
  select(scientific_name, rowS) %>%
  arrange(desc(rowS))
# confused because this doesn't match the dat I've had before... :(

reg_aoo_ind<-reg_aoo %>% filter(area.sqkm!=0) %>%
  group_by(area.type) %>%
  mutate(aoo_ind=(area.sqkm - min(area.sqkm)) / (max(area.sqkm) - min(area.sqkm))) %>%
  arrange(area.type)

write.csv(reg_aoo_ind, '/home/tracidubose/rcs_results/regional_species_shapefiles/Regional_AOO.csv')

reg_aoo_ind<- read.csv('/home/tracidubose/rcs_results/regional_species_shapefiles/Regional_AOO.csv')

# Calculate Climate Sensitivity -----
PATH_buff_climate_val<-'/home/tracidubose/rcs_results/regional_species_shapefiles/Climate_CSV/'
PATH_regional_spp_shp<-'/home/tracidubose/rcs_results/regional_species_shapefiles/'

unique(sub('\\..*','', 
           list.files(paste0(PATH_regional_spp_shp, sub('_sf', '', "Midcontinent"), '/')))) 
unique(sub('\\_.*','',list.files(paste0(PATH_buff_climate_val))))


for(y in usgs_reg_sf[1]){
  for(u in unique(sub('\\..*','', 
                      list.files(paste0(PATH_regional_spp_shp, sub('_sf', '', y), '/'))))){
  extract_mean_env_buff(env.raster = prism_ppt, variable='ppt', 
                        spp=u, native_range = T,
                        buff_file_folder = paste0(PATH_regional_spp_shp, sub('_sf', '', y), '/'), 
                        output_path = paste0(PATH_buff_climate_val, y, '_'))
    extract_mean_env_buff(env.raster = prism_tmin, variable='tmin', 
                          spp=u, native_range = T,
                          buff_file_folder = paste0(PATH_regional_spp_shp, sub('_sf', '', y), '/'), 
                          output_path = paste0(PATH_buff_climate_val, y, '_'))
    extract_mean_env_buff(env.raster = prism_tmax, variable='tmax', 
                          spp=u, native_range = T,
                          buff_file_folder = paste0(PATH_regional_spp_shp, sub('_sf', '', y), '/'), 
                          output_path = paste0(PATH_buff_climate_val, y, '_'))
  }
  print(y)
}

clim_values<-NULL
for(b in list.files(PATH_buff_climate_val,
                    full.name=T)){
  ms<-read.csv(b) %>%
    mutate(region=sub('\\_.*','', substr(b, 72, 85)))
  clim_values<-bind_rows(clim_values, ms)
}
reg_cs_vals<-clim_values %>% group_by(species, region) %>%
  pivot_longer(cols=c('ppt', 'tmax', 'tmin')) %>%
  filter(!is.na(value)) %>%
  group_by(species, region, name) %>%
  dplyr::summarize(mean_value=wtd.mean(value, weight, na.rm = TRUE, normwt = TRUE),
                   sd_value=sqrt(wtd.var(value, weight, na.rm = TRUE, normwt = TRUE))) %>%
  group_by(region, name) %>%
  mutate(var_CS_index = (sd_value-min(sd_value, na.rm=T))/ #calculating the index
           (max(sd_value, na.rm=T)-min(sd_value, na.rm=T))) 

reg_cs_vals %>% filter(is.na(sd_value)) #investigate these - likely just have one point
#Pseudacris maculata
list.files(paste0(PATH_regional_spp_shp, "/Northeast"))
ttest<-read_sf(paste0(PATH_regional_spp_shp, "/Northeast"), layer='Pseudacris maculata')
ggplot()+geom_sf(data=ttest)

regional_rcs<-reg_cs_vals %>% 
  group_by(species, region)%>% 
  dplyr::summarize(CS_val=mean(var_CS_index)) %>%  
  left_join(reg_aoo_ind %>% mutate(region=sub('\\_.*', '', area.type)), 
            by=c('species'='scientific_name', 'region'))%>%
  mutate(RCS=(CS_val+aoo_ind)/2) 

chloropleth_vals<- reg_cs_vals %>% 
  left_join(reg_aoo_ind %>% mutate(region=sub('\\_.*', '', area.type)), 
            by=c('species'='scientific_name', 'region'))%>%
  group_by(region, name) %>%
  dplyr::summarise(n.sp=length(unique(species)),
                   med.area=median(area.sqkm, na.rm=T),
                   med.clim.sd=median(sd_value, na.rm=T)) %>%
  pivot_wider(values_from = med.clim.sd) %>%
  pivot_longer(-region) %>%
  group_by(name) %>%
  mutate(val_scale=scale(value))

states %>% 
  left_join(chloropleth_vals, by=c('usgs_region'='region')) %>%
  ggplot()+
  geom_sf(aes(fill=val_scale))+
  scale_fill_distiller('Scaled (SD) Value', type='div',
                       palette = 'PRGn')+
  facet_wrap(~name)+
  theme_bw()+
  theme(legend.position='bottom')
ggsave('/home/tracidubose/rcs_results/regional_maps_sdval.jpg',
       width=7, height=6)

ntax_region<-reg_aoo_ind %>% group_by(area.type) %>%
  tally() %>%
  mutate(region=sub('\\_.*', '', area.type))

reg_cs_vals %>% 
  bind_rows(reg_aoo_ind %>% 
              mutate(region=sub('\\_.*', '', area.type),
                     name='area.sqkm',
                     log.area=log10(area.sqkm)) %>%
              rename(species='scientific_name', 
                     sd_value='log.area'))%>%
  group_by(region, name) %>% 
  ggplot()+
  geom_density(aes(x=sd_value, fill=region),
               alpha=0.3)+
  scale_fill_viridis_d()+
  facet_wrap(~name, scales='free')
  
taxa_all_regions<-reg_cs_vals %>% 
  group_by(species) %>% 
  dplyr::summarize(n=n()/3) %>%
  arrange(desc(n)) %>% filter(n==6) %>% pull(species)
medRS<-regional_rcs %>% 
  group_by(region) %>% mutate(medRCS=median(RCS, na.rm=T))
regional_rcs %>% 
  group_by(region) %>%
  ggplot()+
  geom_density(aes(x=RCS)) +
  geom_vline(data=medRS, aes(xintercept=medRCS), color='lightgrey')+
  geom_text(data=ntax_region, 
            aes(x=.3, y=.15, label=n))+
  geom_text(data=ntax_region, 
            aes(x=.12, y=.15, label=' spp'))+
  facet_wrap(~region, scales='free')+
  theme_cowplot()+ theme(axis.text.x = element_text(size=9)) + 
  scale_x_reverse(labels=function(x){1-x}, limits=c(1,0))
ggsave('/home/tracidubose/rcs_results/regional_rcs_dist.jpg',
       width=6, height = 4)  

library(ggrepel)
region_shifts<-NULL
reg_cs_vals %>% 
  filter(species %in% taxa_all_regions) %>%
  group_by(species, region)%>% 
  dplyr::summarize(CS_val=mean(var_CS_index)) %>%  
  left_join(reg_aoo_ind %>% mutate(region=sub('\\_.*', '', area.type)), 
            by=c('species'='scientific_name', 'region'))%>%
  mutate(RCS=(CS_val+aoo_ind)/2) %>%
  ggplot()+
  geom_line(aes(x=region, y=RCS, 
                group=species, color=species),
            alpha=0.4)+
  geom_point(aes(x=region, y=RCS), alpha=0.5)+
  #facet_wrap(~name, scales = 'free', ncol=1)+
  scale_color_viridis_d(guide=F)+
  scale_y_reverse(label=function(x){1-x})+
  scale_x_discrete(expand=c(.05,.05))+
  theme_bw()


NW_rcs<-reg_cs_vals %>% 
  #filter(species == 'Anaxyrus woodhousii') %>%
  group_by(species, region)%>% 
  dplyr::summarize(CS_val=mean(var_CS_index)) %>%  
  left_join(reg_aoo_ind) %>%
  mutate(RCS=(CS_val+AOO_index)/2) %>%
  filter(region == 'Northwest_sf') %>%
  left_join(entire_range_rcs) %>%
  select(species, region, RCS, RCS_buff) %>%
  pivot_longer(cols=c(RCS, RCS_buff))%>%
  ggplot()+
  geom_line(aes(x=name, y=value, group=species), color='lightgrey')+
  geom_point(aes(x=name, y=value, color=species), size=2, alpha=0.5)+
  geom_text_repel(data=. %>% filter(name=="RCS"),
                  aes(x='RCS', y=value, label=species),
                  size=3,
                  point.padding=0.1)+
  scale_color_viridis_d(guide=F)+
  scale_y_reverse('', label=function(x){1-x}, lim=c(1,0))+
  scale_x_discrete('', labels=c('Northwest', "Entire Range"))+
  theme_bw()+
  theme(axis.title.y=element_text(size=-3))
library(cowplot)
plot_grid(region_shifts, NW_rcs, rel_widths=c(.5, .3))
ggsave('/home/tracidubose/regional_difs_anuranRCS.jpg',
       width=7, height=4)
