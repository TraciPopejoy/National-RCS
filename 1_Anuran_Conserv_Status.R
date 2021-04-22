# Recording Anuran Conservation Status from three geopolitical units- 'Anuran_Conserv_Status.R'
# Revised by Traci DuBose, 06April2021

library(tidyverse)

#paths for data files
taxa_table<-read.csv('gbif_anuran_taxa_info.csv')
PATH_NatSGCN<-'data/National_SpGCC.csv'
PATH_FedESA<-'data/FWS_Species_Data_Explorer.csv'
counts<-read.csv('rcs_results/occurrence_counts_used.csv') %>% select(-X)

#isolate anuran genus for quick pull from Nat SGCN
anuran_genus<-taxa_table$genus

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
  dplyr::select(scientific_name, nStates_2015, States.2015) %>%
  dplyr::rename(sp_id_gcc='scientific_name') %>%
  #isolate subsp nomen to remove it for easy joining later
  #will produce a lot of NAs and a warning (since most species don't have a subspecies)
  separate(sp_id_gcc, into=c('genus','species','subsp'), sep=" ") %>% 
  rowwise() %>%
  mutate(scientific_name=paste(genus, species)) %>% ungroup()

NatGCC_subset<-NatGCC_list %>% filter(genus %in% c(taxa_table$genus, 'Hyla')) %>% 
  mutate(genus=recode(genus, "Hyla"="Dryophytes"),
         species=recode(species, 
                        'cinerea'='cinereus', 'gratiosa'='gratiosus', 'squirella'='squirellus'),
         taxa=paste(genus, species))

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
         value %in% c(anuran_genus, 'Hyla','Bufo')) %>% #only keep anuran genus 
  mutate(sname=recode(sname, 'Rana chiricahuensis'='Lithobates chiricahuensis',
                      'Bufo houstonensis'='Anaxyrus houstonensis',
                      'Bufo hemiophrys baxteri'= 'Anaxyrus baxteri'))

# Species Status based on IUCN redlist
#install.packages('rredlist')
library(rredlist)

#need an API to use these -- they are free but take days to get
#rredlist::rl_use_iucn() 
#usethis::edit_r_environ()

IUCN_anuran<-data.frame(scientific_name=as.character(),
                        IUCNcategory=as.character(),
                        assess_date=as.character(),
                        population_status=as.character(),
                        stringsAsFactors = F)

for(j in 1:nrow(RCS_index)){
  taxa<-as.character(gsub('\\.',' ',RCS_index$species)[j])
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
IUCN_anuran %>% filter(IUCNcategory=="EW")


con_info<-full_join(IUCN_anuran, esa_anuran_spp, by=c('scientific_name'='sname')) %>%
  full_join(NatGCC_subset, by=c('scientific_name'='taxa')) %>% 
  select(scientific_name, IUCNcategory, ESA.Listing.Status, nStates_2015) %>% 
  mutate(NatSGCN=case_when(nStates_2015 != 0~'in ≥1 SWAP',
                           T~'Not in SWAP')) %>%
  filter(scientific_name %in% counts$final.taxa) %>%
  filter(!duplicated(scientific_name)) %>%
  #replacing NAs from ESA to Not on ESA
  replace_na(list(ESA.Listing.Status = 'Not Listed')) 
write.csv(con_info, 'data/anuran_conservation_info.csv', row.names = F)


# Appendix 1 Table 1 -----
taxa_table %>%
  select(family, species, scientificName, usageKey) %>%
  left_join(counts, by=c('species'='final.taxa')) %>%
  left_join(con_info, by=c('species'='scientific_name')) %>%
  select(-nStates_2015) %>%
  write.csv('rcs_results/Appendix_Table.csv')

# Table 1: N in each group ----
con_info %>% 
  filter(!(scientific_name %in% c('Anaxyrus baxteri', 'Anaxyrus houstonensis'))) %>% 
  pivot_longer(-scientific_name) %>%
  group_by(name, value) %>% tally()



