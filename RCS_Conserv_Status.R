# Investigating how the Anuran RCS aligns with conservation status - 'RCS_Conserv_Status.R'
# Revised by Traci DuBose, 20Sept2020

library(tidyverse)

#paths for data files
PATH_RCS_index<-'results/RCS_table_20200920.csv'
PATH_NatSGCN<-'data/National_SpGCC.csv'
PATH_FedESA<-'data/FWS_Species_Data_Explorer.csv'

#load the RCS index results
RCS_index<-read.csv(PATH_RCS_index)

#isolate anuran genus for quick pull from Nat SGCN
anuran_genus<-unique(gsub(";.*","", gsub("\\ ", ";", RCS_index$scientific_name)))

# Load the lists of species of conservation concern
# National Species of Greatest Conservation Need - https://www1.usgs.gov/csas/swap/
# downloaded 05Aug2020
NatGCC_list<-read.delim(PATH_NatSGCN, 
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
#so problem --- need to make sure all these match the final.taxa/are up to date

# Species on the Endangered Species List - https://ecos.fws.gov/ecp/services
# downloaded 015Aug2020
esa_anuran_spp <- read.csv(PATH_FedESA, stringsAsFactors = F) %>%
  filter(Common.Name!="Lemurs") %>% #long name that was hard to parse
  # scientific name column has strange taxonomic indicators for synonyms
  mutate(sname=Scientific.Name) %>%
  # isolating each word within Scientific.Name
  separate(Scientific.Name, 
           into=c('n1','n2','n3', 'n4', 'n5', 'n6', 'n7','n8'), 
           sep=" ") %>%
  # create long dataframe
  pivot_longer(cols=starts_with('n')) %>%
  filter(!is.na(value),
         name=='n1', #should be the genus for each taxa
         value %in% anuran_genus) #only keep anuran genus 

# Species Status based on IUCN redlist
install.packages('rredlist')
library(rredlist)

#need an API to use these -- they are free but take days to get
rredlist::rl_use_iucn() 
usethis::edit_r_environ()

library(rredlist);library(tidyverse)
IUCN_anuran<-data.frame(scientific_name=as.character(),
                        IUCNcategory=as.character(),
                        assess_date=as.character(),
                        population_status=as.character(),
                        stringsAsFactors = F)

for(j in 1:nrow(RCS_index)){
  taxa<-as.character(gsub('\\.',' ',RCS_index$scientific_name)[j])
  iucn_res<-rl_search(taxa)[[2]]
  IUCN_anuran[j,]$scientific_name<-taxa
  if(is.null(names(iucn_res))){IUCN_anuran[j,2:4]<-'UNK'
  }else{
  IUCN_anuran[j,]$IUCNcategory<-iucn_res$category
  IUCN_anuran[j,]$assess_date<-iucn_res$assessment_date
  IUCN_anuran[j,]$population_status<-iucn_res$population_trend
  }
}
  
#View(IUCN_anuran)
IUCN_anuran %>% filter(IUCNcategory=='UNK')
# Acris blanchardi was recently split from Acris crepitans
# Lithobates kauffeldi was recently moved from Rana
# I can't find either in the IUCN red list.

#for later recoding of the variables
iucn_cat<-c(UNK='Not listed', LC='Least Concern', NT='Near Threatened', VU='Vulnerable', 
  EN='Endangered', CR='Critically Endangered', EW='Extinct in the wild', EX ='Extinct')

### RCS & conservation graph ####
RCS_index_con<-RCS_index %>%
  #only want these two columns
  dplyr::select(scientific_name, RCS_WS) %>%
  #cleaning up the scientific name for easy joining
  mutate(scientific_name=gsub('\\.',' ', scientific_name)) %>%
  #adding in the ESA list
  left_join(esa_anuran_spp, by=c('scientific_name'='sname')) %>%
  #reordering columns
  dplyr::select(scientific_name, RCS_WS, ESA.Listing.Status, Entity.Description, ESA.Listing.Date) %>%
  #adding in the National SGCN list
  left_join(NatGCC_list, by="scientific_name") %>%
  #grouping these into on the list and not on the list
  mutate(NatSGCN=case_when(!is.na(States.2015)~'in >1 SWAP',
                           T~'Not in SWAP')) %>%
  #joining in IUCN data
  left_join(IUCN_anuran, by='scientific_name') %>%
  #replacing NAs from ESA to Not on ESA
  replace_na(list(ESA.Listing.Status = 'Not on ESA')) %>%
  #reordering factors so they plot in an intuitive order (less concern on left)
  mutate(`Endangered Species Act`=factor(ESA.Listing.Status, levels = c('Not on ESA','Threatened','Endangered')),
         `Species of Greatest\nConserv. Need`=factor(NatSGCN, levels=c('Not in SWAP', 'in >1 SWAP')),
         `IUCN Red List`=recode_factor(IUCNcategory, !!!iucn_cat)) %>%
  #keeping only the columns I plan to plot
  dplyr::select(scientific_name, RCS_WS, `Endangered Species Act`,
                `Species of Greatest\nConserv. Need`, `IUCN Red List`) %>%
  #adding a column for species with range issues 
  mutate(RangeIssue=case_when(scientific_name %in% c('Smilisca baudinii','Dryophytes eximius',
                                                    'Leptodactylus fragilis','Smilisca fodiens',
                                                    'Rhinophrynus dorsalis','Anaxyrus hemiophrys',
                                                    'Eleutherodactylus guttilatus')~'range issue',
                             T~'no range issue'))
str(RCS_index_con)

plot_common_elements<-ggplot(data=RCS_index_con)+
  scale_color_manual(guide=F, values=c('grey', 'black'))+
  scale_y_continuous("RCS at HUC12 grain")+
  scale_x_discrete("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=30, hjust=.9))
plot_iucn<-plot_common_elements+
  geom_jitter(aes(x=`IUCN Red List`, y=RCS_WS, color=RangeIssue))+
  geom_boxplot(aes(x=`IUCN Red List`, y=RCS_WS), alpha=0.3,
               outlier.alpha = 0)+
  ggtitle('\nIUCN Red List')
plot_esa<-plot_common_elements+
  geom_jitter(aes(x=`Endangered Species Act`, y=RCS_WS, color=RangeIssue))+
  geom_boxplot(aes(x=`Endangered Species Act`, y=RCS_WS), alpha=0.3,
               outlier.alpha = 0)+
  ggtitle('Endangered\nSpecies Act')
plot_sgcn<-plot_common_elements+
  geom_jitter(aes(x=`Species of Greatest\nConserv. Need`, y=RCS_WS, color=RangeIssue))+
  geom_boxplot(aes(x=`Species of Greatest\nConserv. Need`, y=RCS_WS), alpha=0.3,
               outlier.alpha = 0)+
  ggtitle('Species of Greatest\nConserv. Need')

#plot all in one row
install.packages('cowplot')
library(cowplot)
plot_grid(plot_esa, plot_sgcn, plot_iucn, nrow=1, 
          rel_widths = c(.8,.8,1))
ggsave(paste0('results/RCS_con_', format(Sys.Date(),'%Y%m%d'),
              '.jpg'), width=6, height=4)

# RCS & conservation status quantification ----

# isolating sub species identified in SWAP or ESA ----
subsp.natgcc <- NatGCC_list %>%
  filter(!is.na(subsp),
         genus %in% anuran_genus) %>%
  mutate(sspp=paste(genus, species, subsp)) %>%
  pull(sspp)
subsp.esa <- esa_anuran_spp %>% 
  filter(Entity.Description != 'Wherever found')

occurance.data<-read.csv('data/anuran_occ_all20200708.csv') %>%
  mutate(subspp=case_when(is.na(Taxon)~paste(species, infraspecificEpithet),
                          T~Taxon))
subspp_occ <- occurance.data %>%
  filter(subspp %in% subsp.natgcc)
unique(subspp_occ$subspp)
