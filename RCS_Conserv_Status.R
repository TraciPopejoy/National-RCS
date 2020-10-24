# Investigating how the Anuran RCS aligns with conservation status - 'RCS_Conserv_Status.R'
# Revised by Traci DuBose, 20Sept2020

library(tidyverse)

#paths for data files
PATH_RCS_index<-'rcs_results/RCS_table_20201024.csv'
PATH_RCS_L48_index <-'rcs_results/RCS_table_L48_20201024.csv'
PATH_NatSGCN<-'data/National_SpGCC.csv'
PATH_FedESA<-'data/FWS_Species_Data_Explorer.csv'

#load the RCS index results
RCS_index<-read.csv(PATH_RCS_index)

#isolate anuran genus for quick pull from Nat SGCN
anuran_genus<-unique(gsub(";.*","", gsub("\\ ", ";", RCS_index$species)))

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
         value %in% anuran_genus) %>% #only keep anuran genus 
  mutate(sname=recode(sname, 'Rana chiricahuensis'='Lithobates chiricahuensis'))

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

#for later recoding of the variables
iucn_cat<-c(UNK='Not listed', LC='Least Concern', NT='Near Threatened', VU='Vulnerable', 
  EN='Endangered', CR='Critically Endangered', EW='Extinct in the wild', EX ='Extinct')

### RCS & conservation graph dataframe ####
RCS_index_con<-RCS_Data %>%
  bind_rows(RCS_Data_L48) %>% 
  dplyr::select(species, spatial_extent, RCS_WS, RCS_buff) %>% #only want these columns
  #adding in the ESA list
  left_join(esa_anuran_spp, by=c('species'='sname')) %>%
  #reordering columns
  dplyr::select(species, spatial_extent, RCS_WS, RCS_buff, ESA.Listing.Status, Entity.Description, ESA.Listing.Date) %>%
  #adding in the National SGCN list
  left_join(NatGCC_list, by=c('species'="scientific_name")) %>%
  #grouping these into on the list and not on the list
  mutate(NatSGCN=case_when(!is.na(States.2015)~'in >1 SWAP',
                           T~'Not in SWAP')) %>%
  #joining in IUCN data
  left_join(IUCN_anuran, by=c('species'="scientific_name")) %>%
  #replacing NAs from ESA to Not on ESA
  replace_na(list(ESA.Listing.Status = 'Not on ESA')) %>%
  #reordering factors so they plot in an intuitive order (less concern on left)
  mutate(`Endangered Species Act`=factor(ESA.Listing.Status, levels = c('Not on ESA','Threatened','Endangered')),
         `Species of Greatest\nConserv. Need`=factor(NatSGCN, levels=c('Not in SWAP', 'in >1 SWAP')),
         `IUCN Red List`=recode_factor(IUCNcategory, !!!iucn_cat)) %>%
  #keeping only the columns I plan to plot
  dplyr::select(species, spatial_extent, RCS_WS, RCS_buff, `Endangered Species Act`,
                `Species of Greatest\nConserv. Need`, `IUCN Red List`) %>%
  filter(!duplicated(.))
str(RCS_index_con)

RCS_index_con %>%
  pivot_longer(cols=c('RCS_WS','RCS_buff'), names_to='grain size', values_to='RCS') %>%
  pivot_longer(cols=c(`Endangered Species Act`,
                      `Species of Greatest\nConserv. Need`,
                      `IUCN Red List`), names_to='Geopolitical Group', values_to='Status') %>%
  ggplot()+
  geom_boxplot(aes(x=Status, y=RCS, fill=`Geopolitical Group`))+
  facet_wrap(~spatial_extent*`grain size`, 
             scales = 'free')+
  theme_cowplot()+
  theme(axis.text.x = element_text(size=8, angle=30, hjust=.9),
        legend.position = 'none',)

# RCS & conservation status quantification ----
# differences in mean between conservation status
RCS_index_con_icn<-RCS_index_con %>% #removing small groups
  dplyr::select(-`Endangered Species Act`, -`Species of Greatest\nConserv. Need`) %>%
  filter(!(`IUCN Red List` %in% c('Not listed', 'Extinct in the wild', 
                                  'Extinct'))) %>%
  mutate(`IUCN Red List`=factor(`IUCN Red List`))
# buffered points
b_swap_k<-kruskal.test(formula=RCS_buff~`Species of Greatest\nConserv. Need`,
                       data=RCS_index_con[RCS_index_con$spatial_extent=="entire native NA range",])
b_esa_k<-kruskal.test(formula=RCS_buff~`Endangered Species Act`,
                      data=RCS_index_con[RCS_index_con$spatial_extent=="entire native NA range",])
boxplot(formula=RCS_buff~`IUCN Red List`, data=RCS_index_con_icn)
b_iucn_k<-kruskal.test(formula=RCS_buff~`IUCN Red List`,
                       data=RCS_index_con_icn[RCS_index_con_icn$spatial_extent=="entire native NA range",])

# watershed
ws_swap_k<-kruskal.test(formula=RCS_WS~`Species of Greatest\nConserv. Need`,
                        data=RCS_index_con[RCS_index_con$spatial_extent=="entire native NA range",])
ws_esa_k<-kruskal.test(formula=RCS_WS~`Endangered Species Act`,
                       data=RCS_index_con[RCS_index_con$spatial_extent=="entire native NA range",])
boxplot(formula=RCS_WS~`IUCN Red List`, data=RCS_index_con_icn)
ws_iucn_k<-kruskal.test(formula=RCS_WS~`IUCN Red List`,
                        data=RCS_index_con_icn[RCS_index_con_icn$spatial_extent=="entire native NA range",])

# buffered points
b_swap_k_l48<-kruskal.test(formula=RCS_buff~`Species of Greatest\nConserv. Need`,
                       data=RCS_index_con[RCS_index_con$spatial_extent=="native L48 range",])
b_esa_k_l48<-kruskal.test(formula=RCS_buff~`Endangered Species Act`,
                      data=RCS_index_con[RCS_index_con$spatial_extent=="native L48 range",])
b_iucn_k_l48<-kruskal.test(formula=RCS_buff~`IUCN Red List`,
                       data=RCS_index_con_icn[RCS_index_con_icn$spatial_extent=="native L48 range",])

# watershed
ws_swap_k_l48<-kruskal.test(formula=RCS_WS~`Species of Greatest\nConserv. Need`,
                        data=RCS_index_con[RCS_index_con$spatial_extent=="native L48 range",])
ws_esa_k_l48<-kruskal.test(formula=RCS_WS~`Endangered Species Act`,
                       data=RCS_index_con[RCS_index_con$spatial_extent=="native L48 range",])
ws_iucn_k_l48<-kruskal.test(formula=RCS_WS~`IUCN Red List`,
                        data=RCS_index_con_icn[RCS_index_con_icn$spatial_extent=="native L48 range",])

library(betareg)
summary(betareg(RCS_WS~`IUCN Red List`+spatial_extent, data=RCS_index_con_icn))

# making a table of kruskal wallis outputs
kruskal_res<-sapply(list(ws_esa_k, ws_swap_k, ws_iucn_k,
                         b_esa_k, b_swap_k, b_iucn_k,
                         ws_esa_k_l48, ws_swap_k_l48, ws_iucn_k_l48,
                         b_esa_k_l48, b_swap_k_l48, b_iucn_k_l48), unlist) %>% 
  t() %>% as_tibble() %>% 
  mutate(grain.size=ifelse(grepl('WS', data.name), 'watershed','buffer'),
         spatial_extent=c(rep('entire range', 6), rep('contig US', 6)),
         Chi2=round(as.numeric(`statistic.Kruskal-Wallis chi-squared`),3),
         p=round(as.numeric(p.value), 4),
         star=ifelse(p<0.1, '*', NA)) %>%
  rename(df='parameter.df') %>%
  dplyr::select(spatial_extent, grain.size, data.name, df, Chi2, p, star)
kruskal_res %>% arrange(star)

kruskal_res %>% group_by(data.name, grain.size) %>%
  pivot_wider(names_from = 'spatial_extent', values_from='p') %>%
  summarize(`entire range`=mean(`entire range`, na.rm=T),
            `contig US`=mean(`contig US`, na.rm=T)) %>%
  ggplot()+
  geom_linerange(aes(x=data.name, ymin=`entire range`, ymax=`contig US`))+
  geom_point(aes(x=data.name, y=`entire range`, color='Entire'), size=3)+
  geom_point(aes(x=data.name, y=`contig US`, color='L48'), size=3)+
  scale_x_discrete('')+
  scale_y_continuous('p value', trans='log10')+
  facet_wrap(~grain.size, scales='free',ncol=1,
             strip.position = 'left')+
  coord_flip()+
  theme_cowplot()+
  theme(legend.position='bottom',
        strip.placement = 'outside',
        strip.background = element_rect(fill=NA))

#plot for mean differences
plot_common_elements<-ggplot(data=RCS_index_con)+
  scale_x_discrete("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=30, hjust=.9),
        title = element_text(size=9))
plot_iucn_ws<-plot_common_elements+
  geom_jitter(aes(x=`IUCN Red List`, y=RCS_WS), color='grey', alpha=0.5)+
  geom_boxplot(aes(x=`IUCN Red List`, y=RCS_WS), alpha=0.5,
               outlier.alpha = 0)+
  geom_text(x='Endangered', y=.25,
            label="paste('H'[3]*' = 8.64, p = 0.04')",
            parse=T, size=2)+
  scale_y_continuous('')+
  ggtitle('\nIUCN Red List')+
  theme(axis.title.y = element_text(size=-3))
plot_esa_ws<-plot_common_elements+
  geom_jitter(aes(x=`Endangered Species Act`, y=RCS_WS), color='grey', alpha=0.5)+
  geom_boxplot(aes(x=`Endangered Species Act`, y=RCS_WS), alpha=0.3,
               outlier.alpha = 0)+
  geom_text(x='Threatened', y=.25,
            label="paste('H'[2]*' = 0.31, p = 0.86')",
            parse=T, size=2)+
  scale_y_continuous('Watershed RCS')+
  ggtitle('Endangered\nSpecies Act')
plot_sgcn_ws<-plot_common_elements+
  geom_jitter(aes(x=`Species of Greatest\nConserv. Need`, y=RCS_WS), color='grey', alpha=0.5)+
  geom_boxplot(aes(x=`Species of Greatest\nConserv. Need`, y=RCS_WS), alpha=0.3,
               outlier.alpha = 0)+
  geom_text(x='Not in SWAP', y=.25,
            label="paste('H'[1]*' = .30, p = 0.58')",
            parse=T, size=2)+
  ggtitle('Species of Greatest\nConserv. Need')+
  scale_y_continuous('')+
  theme(axis.title.y = element_text(size=-3))
#buff graphs
plot_iucn_b<-plot_common_elements+
  geom_jitter(aes(x=`IUCN Red List`, y=RCS_buff), color='grey', alpha=0.5)+
  geom_boxplot(aes(x=`IUCN Red List`, y=RCS_buff), alpha=0.5,
               outlier.alpha = 0)+
  geom_text(x='Endangered', y=.25,
            label="paste('H'[3]*' = 6.23, p = 0.10')",
            parse=T, size=2)+
  scale_y_continuous('')+
  ggtitle('\nIUCN Red List')+
  theme(axis.title.y = element_text(size=-3))
plot_esa_b<-plot_common_elements+
  geom_jitter(aes(x=`Endangered Species Act`, y=RCS_buff), color='grey', alpha=0.5)+
  geom_boxplot(aes(x=`Endangered Species Act`, y=RCS_buff), alpha=0.3,
               outlier.alpha = 0)+
  geom_text(x='Threatened', y=.25,
            label="paste('H'[2]*' = 0.79, p = 0.67')",
            parse=T, size=2)+
  scale_y_continuous('1 km buffer RCS')+
  ggtitle('Endangered\nSpecies Act')
plot_sgcn_b<-plot_common_elements+
  geom_jitter(aes(x=`Species of Greatest\nConserv. Need`, y=RCS_buff), color='grey', alpha=0.5)+
  geom_boxplot(aes(x=`Species of Greatest\nConserv. Need`, y=RCS_buff), alpha=0.3,
               outlier.alpha = 0)+
  geom_text(x='Not in SWAP', y=.25,
            label="paste('H'[1]*' = 0.14, p = 0.71')",
            parse=T, size=2)+
  scale_y_continuous('')+
  ggtitle('Species of Greatest\nConserv. Need')+
  theme(axis.title.y = element_text(size=-3))

#plot all in one row
library(cowplot)
plot_grid(plot_esa_ws, plot_sgcn_ws, plot_iucn_ws,
          plot_esa_b, plot_sgcn_b, plot_iucn_b,
          nrow=2, rel_widths = c(.85,.85,1), labels=c('A','','','B','',''))
ggsave(paste0('rcs_results/figures/RCS_con_means', format(Sys.Date(),'%Y%m%d'),'.jpg'), 
       width=6, height=6)

# RCS to predict whether considered threatened or not ----
#stat notes & examples
# using Practical Guide to Logistic Regression by Joseph M. Hilbe (2015) as a guide
# also https://m-clark.github.io/docs/GAMS.pdf for explanation of GAMS
library(mgcv)

rcs<-RCS_index_con %>%
  mutate(ESA_bin=case_when(`Endangered Species Act`=='Not on ESA'~0,
                           T~1),
         State_bin=case_when(`Species of Greatest\nConserv. Need`=='in >1 SWAP'~1,
                             T~0),
         Int_bin=case_when(`IUCN Red List` %in% c('Not listed', 'Least Concern', 'Near Threatened')~0,
                           T~1))
hist(rcs$RCS_WS)

esagam <- gam(ESA_bin~s(RCS_WS), family=binomial, data=rcs)
summary(esagam);plot(esagam)
summary(esagam)$s.table
esaglm<-glm(ESA_bin~RCS_WS, family=binomial, data=rcs)
summary(esaglm)
AIC(esagam, esaglm)

stagam<-gam(State_bin~s(RCS_WS), family=binomial, data=rcs) 
summary(stagam);plot(stagam) #edf close to 1, ok to do logistic regression on?
summary(stagam)$s.table
staglm<-glm(State_bin~RCS_WS, family=binomial, data=rcs)
summary(staglm)
AIC(stagam, staglm)

intgam<-gam(Int_bin~s(RCS_WS), family=binomial, data=rcs) 
summary(intgam);plot(intgam)
summary(intgam)$s.table
intglm<-glm(Int_bin~RCS_WS, family=binomial, data=rcs)
summary(intglm)
AIC(intgam, intglm)

log.p.vals<-rbind(summary(intgam)$s.table,
          summary(esagam)$s.table,
          summary(stagam)$s.table) %>%
  as_tibble() %>%
  mutate(name=c('IUCN','ESA','SWAP'))

new_data_fake<-data.frame(State_bin=1,
                          ESA_bin=1,
                          Int_bin=1,
                          RCS_WS=seq(0,1,by=0.01))
loking<-new_data_fake %>%
  mutate('IUCN'=predict(intgam, newdata=new_data_fake,type="response"),
         'ESA'=predict(esagam, newdata=new_data_fake,type="response"),
         'SWAP'=predict(stagam, newdata=new_data_fake,type="response")) %>%
  pivot_longer(cols = c('ESA',"SWAP", "IUCN"))

ggplot()+
  geom_line(data=loking, 
            aes(x=RCS_WS, y=value, color=name), size=2)+
  geom_text(aes(x=c(.95, .87, .7), 
                y=c(.065, .4, .82), 
                label=c('ESA','IUCN*','State')))+
  scale_x_continuous('watershed RCS index')+
  scale_y_continuous('Prob. of being listed')+
  scale_color_viridis_d(option='E', end=.8, guide=F)+
  theme_cowplot()+theme(legend.position = 'top')+
  ggtitle('Logistic Regression with cubic spline')

rcs_ddddd<-rcs %>%
  rename('ESA'=ESA_bin, 'SWAP'=State_bin, 'IUCN'=Int_bin) %>%
  pivot_longer(cols = c('ESA',"SWAP", "IUCN"))
ggplot()+
  geom_line(data=loking, 
            aes(x=RCS_WS, y=value, color=name), size=1.5)+
  geom_point(data=rcs_ddddd, aes(x=RCS_WS, y=value, 
                                 fill=as.factor(value)),
             alpha=0.3, position=position_jitter(height=.02),
             shape=21)+
  geom_text(data=log.p.vals, x=0, y=1,
            aes(label=paste('p = ', round(`p-value`,2))),
            hjust=0)+
  scale_x_continuous('watershed RCS')+
  scale_y_continuous('Prob. of being listed')+
  scale_colour_viridis_d('List',end=.8, option='E', guide=F)+
  scale_fill_manual('Listed as\nThreatened',values=c('white', 'black'),
                    labels=c('No','Yes'))+
  facet_wrap(~name, ncol=1)+
  theme_cowplot()+  theme(legend.position = 'bottom',
                          axis.text.x = element_text(size=10),
                          axis.text.y = element_text(size=10),
                          legend.box="vertical",
                          legend.justification = 'center') +
  guides(fill=guide_legend(nrow=2,
                           override.aes=list(alpha=1)))
ggsave('rcs_results/figures/logreg_output.jpg', width=3.5, height=7)
summary(esagam)
summary(intgam)
summary(stagam)



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


### which species are unlisted but high RCS ----
esa_spp_nc<-RCS_index_con %>%
  filter(`Endangered Species Act`=='Not on ESA') %>%
  arrange(desc(RCS_WS)) %>% slice(1:7) %>%
  dplyr::select(species, RCS_WS, RCS_buff)
swap_spp_nc<-RCS_index_con %>%
  filter(`Species of Greatest\nConserv. Need`=='Not in SWAP') %>%
  arrange(desc(RCS_WS)) %>% slice(1:7) %>%
  dplyr::select(species, RCS_WS, RCS_buff)
iucn_spp_nc<-RCS_index_con %>%
  filter(!(`IUCN Red List` %in% c('Not listed', 'Least Concern', 'Near Threatened'))) %>%
  arrange(desc(RCS_WS)) %>% slice(1:7) %>%
  dplyr::select(species, RCS_WS, RCS_buff)

spp_to_research<-RCS_index_con %>% 
  filter(species %in% (bind_rows(esa_spp_nc, swap_spp_nc, iucn_spp_nc) %>%
                         filter(!duplicated(.)) %>% pull(species))) %>%
  arrange(desc(species))
write.csv(spp_to_research, 'rcs_results/spp_to_research.csv')
  