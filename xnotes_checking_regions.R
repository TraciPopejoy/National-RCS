#looking at regions
library(sf)
regions<-NULL
for(u in unique(sub('\\..*','', list.files('Region Shapefiles/')))){
  beep<-read_sf('Region Shapefiles', u)
  regions<-bind_rows(regions,beep)
}
regions
ggplot()+geom_sf(data = regions, aes(fill=Name), alpha=0.5)

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
                               ID %in% RockyMountain ~ 'Rocky Mountain'))
ggplot()+geom_sf(data=states, aes(fill=usgs_region))

NatGCC_list<-read.delim('data/National_SpGCC.csv', 
                        sep="|", header=T, stringsAsFactors = F) %>%
  filter(Taxonomic.Group=="Amphibians") %>%
  dplyr::select(-Taxonomic.Group) %>% 
  #renaming columns for ease of use
  rename(scientific_name='Scientific.Name', nStates_2005='X..of.States.2005',
         nStates_2015='X..of.States.2015') %>%
  filter(nStates_2015>0) %>% #removing species that were not found on any state wildlife plan
  #only keeping the columns 
  dplyr::select(scientific_name, States.2015) %>%
  dplyr::rename(sp_id_gcc='scientific_name') %>%
  #isolate subsp nomen to remove it for easy joining later
  #will produce a lot of NAs and a warning (since most species don't have a subspecies)
  separate(sp_id_gcc, into=c('genus','species','subsp'), sep=" ") %>% 
  rowwise() %>%
  mutate(scientific_name=paste(genus, species)) 
anura_genera <- RCS_index %>%
  select(species) %>%
  separate(species, c('genus', 'species')) %>%
  pull(genus) %>% unique() %>%
  c("Hyla", "Bufo")

SGCN_usreg <- NatGCC_list %>% 
  ungroup() %>%
  filter(genus %in% anura_genera,
         scientific_name %in% RCS_index$species,
         !duplicated(.))%>%
  separate(States.2015, into=paste0('x', 1:19), sep=',') %>%
  pivot_longer(cols=starts_with('x')) %>%
  filter(!is.na(value),
         value != "Alaska") %>%
  mutate(list_state=tolower(value),
         usgs_region=case_when(list_state %in% Northwest ~ 'Northwest',
                               list_state %in% Northeast ~ 'Northeast',
                               list_state %in% Southwest ~ 'Southwest',
                               list_state %in% Southeast ~ 'Southeast',
                               list_state %in% Midcontintent ~ 'Midcontinent',
                               list_state %in% RockyMountain ~ 'Rocky Mountain')) %>%
  select(usgs_region, scientific_name) %>%
  mutate(presence=1) %>%
  pivot_wider(names_from = 'usgs_region', values_from='presence', values_fn=sum,
              values_fill=0) %>%
  left_join(NatGCC_list %>% select(scientific_name, States.2015))
SGCN_usreg %>% 
  select(-scientific_name, -States.2015) %>%
  summarize_all(sum, na.rm=T)
states %>% 
  mutate(state.area.sqm=st_area(.)) %>% 
  group_by(usgs_region) %>%
  summarize(reg.area.sqkm=sum(state.area.sqm)/1e6) %>%
  arrange(reg.area.sqkm) 

left_join(SGCN_usreg, AOOs) %>%
  group_by(scientific_name) %>%
  summarize(nStates=sum(`Southwest`:`Midcontinent`, na.rm = T))
View(SGCN_usreg)

library(scales); library(ggrepel)
read.delim(PATH_NatSGCN, 
           sep="|", header=T, stringsAsFactors = F) %>%
  filter(Taxonomic.Group=="Amphibians") %>%
  dplyr::select(-Taxonomic.Group) %>% 
  #renaming columns for ease of use
  rename(scientific_name='Scientific.Name', nStates_2005='X..of.States.2005',
         nStates_2015='X..of.States.2015') %>%
  filter(nStates_2015>0) %>%
  select(scientific_name, nStates_2015) %>%
  filter(scientific_name %in% RCS_index$species) %>%
  left_join(RCS_index, by=c('scientific_name'='species')) %>%
  mutate(n_states=as.numeric(paste(nStates_2015))) %>%
  select(scientific_name, nStates_2015, n_states, watershed, RCS_WS) %>%
  ggplot()+
  geom_point(aes(x=n_states, y=watershed, size=RCS_WS),
             color='darkgrey', alpha=0.3)+
  geom_text_repel(data=. %>% ungroup() %>%
                    slice(1:4),
                  aes(x=n_states, y=watershed, label=scientific_name),
                  size=3)+
  scale_y_log10(expression('Area of Occupancy, km '^2), 
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  scale_x_continuous('Times listed on a SWAP',
                     breaks=c(1,5,10,15,20),
                     limits=c(1,20))+
  scale_size('RCS')+
  theme_cowplot()+
  theme(legend.position='bottom',
        legend.justification = 'center')
ggsave('rcs_results/figures/sup_aoo_sgcn.jpg', width=4, height = 4)
