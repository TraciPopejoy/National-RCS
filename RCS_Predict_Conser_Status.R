# Investigating how the Anuran RCS aligns with conservation status - 'RCS_Predict_Conser_Status.R'
# Revised by Traci DuBose, 06April2021

library(tidyverse);library(cowplot)

#paths for data files
PATH_RCS_index<-'rcs_results/RCS_table_20210405.csv'
PATH_RCS_L48_index <-'rcs_results/RCS_table_L48_20210405.csv'

#load the Conservation Information
con_info<-read.csv('data/anuran_conservation_info.csv')
#for later recoding of the variables
iucn_cat<-c(UNK='Not listed', LC='Least Concern',
            NT='Near Threatened', VU='Vulnerable', 
            EN='Endangered', CR='Critically Endangered', 
            EW='Endangered', #changed for A.boreaus
            EX ='Extinct')

### RCS & conservation graph dataframe ####
RCS_index_con<-read.csv(PATH_RCS_index) %>%
  bind_rows(read.csv(PATH_RCS_L48_index)) %>% 
  dplyr::select(species, spatial_extent, RCS_WS, RCS_buff) %>% #only want these columns
  left_join(con_info, by=c('species'='scientific_name')) %>% # join with the conservation information
  #reordering factors so they plot in an intuitive order (less concern on left)
  mutate(`Endangered Species Act`=factor(ESA.Listing.Status, levels = c('Not Listed','Threatened','Endangered')),
         `Species of Greatest Conserv. Need`=factor(NatSGCN, levels=c('Not in SWAP', 'in ≥1 SWAP')),
         `IUCN Red List`=recode_factor(IUCNcategory, !!!iucn_cat)) %>%
  #keeping only the columns I plan to plot
  dplyr::select(species, spatial_extent, RCS_WS, RCS_buff, `Endangered Species Act`,
                `Species of Greatest Conserv. Need`, `IUCN Red List`) %>%
  filter(!duplicated(.))
str(RCS_index_con)

# RCS & conservation status quantification ----
# differences in mean between conservation status
RCS_index_con_icn<-RCS_index_con %>% #removing small groups
  dplyr::select(-`Endangered Species Act`, -`Species of Greatest Conserv. Need`) %>%
  filter(!(`IUCN Red List` %in% c('Not listed', 'Extinct in the wild', 
                                  'Extinct'))) %>%
  mutate(`IUCN Red List`=factor(`IUCN Red List`))
# buffered points
b_swap_k<-kruskal.test(formula=RCS_buff~`Species of Greatest Conserv. Need`,
                       data=RCS_index_con[RCS_index_con$spatial_extent=="entire native NA range",])
b_esa_k<-kruskal.test(formula=RCS_buff~`Endangered Species Act`,
                      data=RCS_index_con[RCS_index_con$spatial_extent=="entire native NA range",])
boxplot(formula=RCS_buff~`IUCN Red List`, data=RCS_index_con_icn)
b_iucn_k<-kruskal.test(formula=RCS_buff~`IUCN Red List`,
                       data=RCS_index_con_icn[RCS_index_con_icn$spatial_extent=="entire native NA range",])

# watershed
ws_swap_k<-kruskal.test(formula=RCS_WS~`Species of Greatest Conserv. Need`,
                        data=RCS_index_con[RCS_index_con$spatial_extent=="entire native NA range",])
ws_esa_k<-kruskal.test(formula=RCS_WS~`Endangered Species Act`,
                       data=RCS_index_con[RCS_index_con$spatial_extent=="entire native NA range",])
boxplot(formula=RCS_WS~`IUCN Red List`, data=RCS_index_con_icn)
ws_iucn_k<-kruskal.test(formula=RCS_WS~`IUCN Red List`,
                        data=RCS_index_con_icn[RCS_index_con_icn$spatial_extent=="entire native NA range",])

# buffered points
b_swap_k_l48<-kruskal.test(formula=RCS_buff~`Species of Greatest Conserv. Need`,
                           data=RCS_index_con[RCS_index_con$spatial_extent=="native L48 range",])
b_esa_k_l48<-kruskal.test(formula=RCS_buff~`Endangered Species Act`,
                          data=RCS_index_con[RCS_index_con$spatial_extent=="native L48 range",])
b_iucn_k_l48<-kruskal.test(formula=RCS_buff~`IUCN Red List`,
                           data=RCS_index_con_icn[RCS_index_con_icn$spatial_extent=="native L48 range",])

# watershed
ws_swap_k_l48<-kruskal.test(formula=RCS_WS~`Species of Greatest Conserv. Need`,
                            data=RCS_index_con[RCS_index_con$spatial_extent=="native L48 range",])
ws_esa_k_l48<-kruskal.test(formula=RCS_WS~`Endangered Species Act`,
                           data=RCS_index_con[RCS_index_con$spatial_extent=="native L48 range",])
ws_iucn_k_l48<-kruskal.test(formula=RCS_WS~`IUCN Red List`,
                            data=RCS_index_con_icn[RCS_index_con_icn$spatial_extent=="native L48 range",])

# Table S2: kruskal wallis outputs ----
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
kruskal_res %>% arrange(star, p)

library(dunn.test)
d_test_data<-RCS_index_con_icn %>%
  filter(spatial_extent=='entire native NA range')
dunn.test(d_test_data$RCS_WS,
          d_test_data$`IUCN Red List`,
          method="bonferroni", alpha=0.1)

#Figure 2: mean differences among conservation status ----
mean_rcs_plot_facet<-RCS_index_con %>%
  dplyr::select(-RCS_buff) %>%
  pivot_longer(cols=c("IUCN Red List", "Endangered Species Act",
                      "Species of Greatest Conserv. Need")) %>%
  mutate(value=factor(value,
                      levels=c("Not in SWAP","Not Listed","Least Concern",
                               "Near Threatened",
                               "in =1 SWAP","Vulnerable","Threatened",
                               "Endangered"),
                      labels=c("Not in SWAP","Not\nListed","Least\nConcern",
                               "Near\nThreatened",
                               "in \u2265 1 SWAP",
                               "Vulnerable","Threatened", "Endangered")),
         name=factor(name, levels=c("IUCN Red List",
                                    "Endangered Species Act",
                                    "Species of Greatest Conserv. Need"))) %>%
  filter(!is.na(value)) %>% 
  ggplot()+
  geom_boxplot(aes(x=value, y=RCS_WS, fill=spatial_extent))+
  facet_wrap(~name, ncol=1, scales='free_x')+
  scale_y_continuous('Watershed RCS')+
  scale_x_discrete("Conservation Status")+
  scale_color_manual(values=c('darkgrey','lightgrey'),
                     aesthetics=c('fill','color'), guide=F)+
  theme_bw()

plot_common_elements<-ggplot(data=RCS_index_con)+
  scale_y_continuous('Watershed RCS')+
  scale_color_manual(values=c('darkgrey','lightgrey'),
                     aesthetics=c('fill','color'), guide=F)+
  theme_bw()+
  theme(title = element_text(size=9))

plot_iucn_ws<-plot_common_elements+
  #geom_jitter(aes(x=`IUCN Red List`, y=RCS_WS, 
  #                color=spatial_extent, group=spatial_extent))+
  geom_boxplot(aes(x=`IUCN Red List`, y=RCS_WS,
                   fill=spatial_extent))+
  scale_x_discrete("",
                   labels=c("Not\nListed","Least\nConcern",
                            "Near\nThreatened","Vulnerable",
                            "Endangered"))+
  ggtitle(' IUCN Red List')
plot_esa_ws<-plot_common_elements+
  #geom_jitter(aes(x=`Endangered Species Act`, y=RCS_WS), color='grey', alpha=0.5)+
  geom_boxplot(aes(x=`Endangered Species Act`, y=RCS_WS,
                   fill=spatial_extent))+
  #geom_text(x='Threatened', y=.25,
  #          label="paste('H'[2]*' = 0.31, p = 0.86')",
  #          parse=T, size=2)+
  scale_x_discrete("",
                   labels=c("Not\nListed","Threatened","Endangered"))+
  ggtitle('Endangered Species Act')
plot_sgcn_ws<-plot_common_elements+
  #geom_jitter(aes(x=`Species of Greatest Conserv. Need`, y=RCS_WS), color='grey', alpha=0.5)+
  geom_boxplot(aes(x=`Species of Greatest Conserv. Need`, 
                   y=RCS_WS,
                   fill=spatial_extent))+
  #geom_text(x='Not in SWAP', y=.25,
  #          label="paste('H'[1]*' = .30, p = 0.58')",
  #          parse=T, size=2)+
  ggtitle('Species of Greatest Conserv. Need')+
  scale_x_discrete("",
                   labels=c("Not in SWAP",
                            expression("in ">=1*" SWAP")))+
  theme(axis.title.y = element_text(size=-3))

#plot all in one row
mean_rcs_plot<-plot_grid(plot_iucn_ws, plot_esa_ws, plot_sgcn_ws, 
                         ncol=1, labels="AUTO")

RCS_index_con %>%
  filter(spatial_extent=='entire native NA range',
         `Endangered Species Act`=='Endangered') %>%
  pull(RCS_WS) %>% quantile(0.7)

RCS_index_con %>%
  filter(spatial_extent=='entire native NA range',
         RCS_WS >= 0.816)

# Table 2: kruskal wallis median results ----
RCS_index_con %>% 
  mutate(spatial_extent=recode(spatial_extent, 
                               'entire native NA range'='entire range',
                               'native L48 range'='contig US')) %>%
  select(species, spatial_extent, RCS_WS, RCS_buff, 
         `Endangered Species Act`, `Species of Greatest Conserv. Need`, `IUCN Red List`)%>%
  rename(watershed='RCS_WS', buffer='RCS_buff') %>%
  pivot_longer(cols=c('watershed','buffer'), names_to='grain.size') %>% 
  pivot_longer(cols=c(`Endangered Species Act`, `Species of Greatest Conserv. Need`, `IUCN Red List`), 
               names_to='con.group', values_to='con.stat') %>% 
  mutate(con.stat=as.character(con.stat),
         con.stat.good=recode(con.stat, 
                              'Endangered'='conserved', 'Not Listed'='not conserved', 'Threatened'='conserved',
                              'Extinct in the wild'='exclude', 'Least Concern'='not conserved', 'Near Threatened'='not conserved',
                              'Not listed'='exclude', 'Vulnerable'='conserved', 'in =1 SWAP'='conserved',
                              'Not in SWAP'='not conserved')) %>%
  group_by(spatial_extent, grain.size, con.group, con.stat.good) %>%
  summarize(med.RCS=round(median(value),3)) %>% 
  pivot_wider(names_from = con.stat.good, values_from = med.RCS) %>%
  ungroup() %>%
  select(-exclude) %>%
  left_join(kruskal_res  %>% 
              mutate(con.group=case_when(substr(data.name, 5,6)=='WS' ~ substr(data.name, 11, nchar(data.name)),
                                         substr(data.name, 5,6)=='bu' ~ substr(data.name, 13, nchar(data.name)))))%>%
  select(-data.name, -star) %>%
  mutate(con.group=factor(con.group, levels=c('IUCN Red List', 'Endangered Species Act', 
                                              'Species of Greatest Conserv. Need'))) %>%
  arrange(desc(spatial_extent), desc(grain.size), con.group) %>%
  write.csv('rcs_results/kruskal_wallace_res_Table.csv')

# RCS to predict whether considered threatened or not ----
# using Practical Guide to Logistic Regression by Joseph M. Hilbe (2015) as a guide
# also https://m-clark.github.io/docs/GAMS.pdf for explanation of GAMS
library(mgcv)

# Split data and make conservation status a binary variable ---
rcs_na<-RCS_index_con %>%
  mutate(ESA_bin=case_when(`Endangered Species Act`=='Not Listed'~0,
                           T~1),
         State_bin=case_when(`Species of Greatest Conserv. Need`=='in =1 SWAP'~1,
                             T~0),
         Int_bin=case_when(`IUCN Red List` %in% c('Not listed', 'Least Concern', 'Near Threatened')~0,
                           T~1),
         ESA_cat=as.numeric(`Endangered Species Act`)) %>%
  filter(spatial_extent=='entire native NA range')
rcs_l48<-RCS_index_con %>%
  mutate(ESA_bin=case_when(`Endangered Species Act`=='Not Listed'~0,
                           T~1),
         State_bin=case_when(`Species of Greatest Conserv. Need`=='in =1 SWAP'~1,
                             T~0),
         Int_bin=case_when(`IUCN Red List` %in% c('Not listed', 'Least Concern', 'Near Threatened')~0,
                           T~1)) %>%
  filter(spatial_extent=='native L48 range')

# Run logistic regression models ---
esagam_er <- gam(ESA_bin~s(RCS_WS), family=binomial, data=rcs_na)
esaglm_er <- glm(ESA_bin~RCS_WS, family=binomial, data=rcs_na)
#esagam_oc_er<-gam(ESA_cat~RCS_WS, family=ocat(R=3),
#                  data=data.frame(rcs_na))
#summary(esagam_oc_er)
summary(esagam_er)
AIC(esagam_er, esaglm_er) #, esagam_oc_er)
gam.check(esagam_er)
esagam_l4 <- gam(ESA_bin~s(RCS_WS), family=binomial, data=rcs_l48)
esaglm_l4 <- glm(ESA_bin~RCS_WS, family=binomial, data=rcs_l48)
AIC(esagam_l4, esaglm_l4)
gam.check(esagam_l4)
deviance(esagam_er); deviance(esagam_l4)

intgam_er <- gam(Int_bin~s(RCS_WS), family=binomial, data=rcs_na)
gam.check(intgam_er)
intgam_l4 <- gam(Int_bin~s(RCS_WS), family=binomial, data=rcs_l48)
gam.check(intgam_l4)
deviance(intgam_er); deviance(intgam_l4)

stagam_er <- gam(State_bin~s(RCS_WS), family=binomial, data=rcs_na)
gam.check(stagam_er)
stagam_l4 <- gam(State_bin~s(RCS_WS), family=binomial, data=rcs_l48)
gam.check(stagam_l4)
deviance(stagam_er); deviance(stagam_l4)

# Table S4 ----
log.p.vals<-rbind(summary(intgam_er)$s.table,
                  summary(esagam_er)$s.table,
                  summary(stagam_er)$s.table,
                  summary(intgam_l4)$s.table,
                  summary(esagam_l4)$s.table,
                  summary(stagam_l4)$s.table) %>%
  as_tibble() %>%
  mutate(name=rep(c('IUCN','ESA','SWAP'),2),
         spat_extent=rep(c('entire', 'US'), each=3))

# Plot logistic regression predictions ----
new_data_fake<-data.frame(State_bin=1,
                          ESA_bin=1,
                          Int_bin=1,
                          RCS_WS=seq(0,1,by=0.01))
predicted_data<-new_data_fake %>%
  mutate('IUCN_er'=predict(intgam_er, newdata=new_data_fake,type="response"),
         'ESA_er'=predict(esagam_er, newdata=new_data_fake,type="response"),
         'SGCN_er'=predict(stagam_er, newdata=new_data_fake,type="response"),
         'IUCN_l4'=predict(intgam_l4, newdata=new_data_fake,type="response"),
         'ESA_l4'=predict(esagam_l4, newdata=new_data_fake,type="response"),
         'SGCN_l4'=predict(stagam_l4, newdata=new_data_fake,type="response"))%>%
  pivot_longer(cols = -c('State_bin',"ESA_bin", "Int_bin","RCS_WS"))%>% 
  rename(model='name') %>%
  mutate(name=sub("_.*","",model),
         spatial_extent=recode(sub(".*_","",model),
                               er='entire native NA range',
                               l4='native L48 range')) %>%
  mutate(geopol=factor(name, levels=c("IUCN", "ESA","SGCN"),
                       labels=c("IUCN Red List",
                                "Endangered Species Act",
                                "Species of Greatest Conserv. Need")))

ggplot()+
  geom_line(data=predicted_data, 
            aes(x=RCS_WS, y=value, color=name,
                linetype=name), size=1.5)+
  geom_text(aes(x=c(1.05, 1.06, 1.05), 
                y=c(.167, .55, .93), 
                label=c('ESA','IUCN*','State')))+
  scale_linetype_manual(values=c(1,1,2,2,4,4))+
  scale_x_continuous('watershed RCS index', 
                     breaks=seq(0,1,.25),expand=c(0.05,.05))+
  scale_y_continuous('Prob. of being listed')+
  scale_color_viridis_d(option='E')+
  theme_cowplot()+theme(legend.position = 'top')+
  ggtitle('Logistic Regression with cubic spline')

tall.log.pred.p<-rcs_na %>% bind_rows(rcs_l48) %>%
  rename('ESA'=ESA_bin, 'SGCN'=State_bin, 'IUCN'=Int_bin) %>%
  pivot_longer(cols = c('ESA',"SGCN", "IUCN")) %>%
  mutate(geopol=factor(name, levels=c("IUCN", "ESA","SGCN"),
                       labels=c("IUCN Red List",
                                "Endangered Species Act",
                                "Species of Greatest Conserv. Need"))) %>%#View()
  ggplot()+
  geom_line(data=predicted_data, 
            aes(x=RCS_WS, y=value, 
                color=spatial_extent), size=1.5)+
  geom_point(aes(x=RCS_WS, y=value, 
                 fill=as.factor(value)),
             alpha=0.3, position=position_jitter(height=.02),
             shape=21)+
  #geom_text(data=log.p.vals, x=0, y=1,
  #          aes(label=paste('p = ', round(`p-value`,2))),
  #          hjust=0)+
  scale_x_continuous('watershed RCS')+
  scale_y_continuous('Probability of being listed')+
  scale_color_manual(values=c('darkgrey','lightgrey'), guide=F)+
  scale_fill_manual('Listed as Threatened',values=c('white', 'black'),
                    labels=c('No','Yes'), guide=F)+
  facet_wrap(~geopol, ncol=1)+
  theme_bw()+  theme(legend.position = 'bottom',
                     axis.text.x = element_text(size=10),
                     axis.text.y = element_text(size=10),
                     legend.justification = 'center') #+
#guides(fill=guide_legend(nrow=2,
#                         override.aes=list(alpha=1)))

plot_grid(mean_rcs_plot_facet, tall.log.pred.p, nrow=1,
          labels=c("C","D"))

ggsave("rcs_results/figures/Fig2_mean_log.jpeg",
       width=6.5, height=7)

# which are likely under served ----
# ESA top 30th percentile
RCS_index_con %>%
  filter(`Endangered Species Act`=="Endangered",
         spatial_extent=='entire native NA range') %>%
  pull(RCS_WS) %>% quantile(.7)
#cut off of 0.865 for 50%
#cut off of 0.967 for 75%
# underserved species based on top 30th percentile
RCS_index_con %>% filter(RCS_WS >= .816,
                         `Endangered Species Act` != 'Endangered',
                         spatial_extent=='entire native NA range')
#if only conserve same %
RCS_index_con %>%
  filter(spatial_extent == 'entire native NA range') %>%
  arrange(desc(RCS_WS)) %>%
  slice(1:5) %>% select(species, RCS_WS, `Endangered Species Act`) 

# IUCN top 30th percentile
RCS_index_con %>%
  filter(`IUCN Red List`=="Endangered",
         spatial_extent=='entire native NA range') %>%
  pull(RCS_WS) %>% quantile(.7) #0.9202
# underserved species based on top 30th percentile
RCS_index_con %>% filter(RCS_WS >= .892,
                         `IUCN Red List` != 'Endangered',
                         spatial_extent=='entire native NA range') %>%
  arrange(desc(RCS_WS))
#if only conserve same %
RCS_index_con %>%
  filter(spatial_extent == 'entire native NA range') %>%
  arrange(desc(RCS_WS)) %>%
  slice(1:7) %>% select(species, RCS_WS, `IUCN Red List`)

# State top 30th percentile
RCS_index_con %>%
  filter(`Species of Greatest Conserv. Need`=='in =1 SWAP',
         spatial_extent=='entire native NA range') %>%
  pull(RCS_WS) %>% quantile(.7) #.8368
# underserved species based on top 30th percentile
RCS_index_con %>% filter(RCS_WS >= .8368,
                         `Species of Greatest Conserv. Need`!='in =1 SWAP',
                         spatial_extent=='entire native NA range') %>%
  arrange(desc(RCS_WS))
#if only conserve same %
RCS_index_con %>%
  filter(spatial_extent == 'entire native NA range') %>%
  arrange(desc(RCS_WS)) %>%
  slice(1:7) %>% select(species, RCS_WS, `IUCN Red List`)




### which species are unlisted but high RCS ----
esa_spp_nc<-RCS_index_con %>%
  filter(`Endangered Species Act`=='Not on ESA') %>%
  arrange(desc(RCS_WS)) %>% slice(1:7) %>%
  dplyr::select(species, RCS_WS, RCS_buff)
swap_spp_nc<-RCS_index_con %>%
  filter(`Species of Greatest Conserv. Need`=='Not in SWAP') %>%
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

# Look at candidate species ----
# source: https://ecos.fws.gov/ecp/report/adhoc-creator?catalogId=species&reportId=species&columns=%2Fspecies@cn,sn,status,desc,listing_date&sort=%2Fspecies@cn%20asc;%2Fspecies@sn%20asc
# found here: https://ecos.fws.gov/ecp/species-reports > Data Explorers > FWS Species
# definition list: https://www.fws.gov/endangered/about/listing-status-codes.html#:~:text=Resolved%20Taxon%20(RT)%20%2D%20Species,published%20in%20the%20Federal%20Register.
test <- read.csv('C:/Users/Owner/Downloads/FWS_Species_Data_Explorer.csv')
head(test)
unique(test$ESA.Listing.Status)

FWS_spp<-test %>% 
  mutate(sname=Scientific.Name) %>%
  separate(Scientific.Name, 
           into=c('n1','n2','n3', 'n4'), 
           sep=" ") %>%
  filter(n1 %in% c(anuran_genus,'Hyla')) %>% arrange(Common.Name) %>%
  anti_join(esa_anuran_spp, by=c('Common.Name', 'sname')) %>%
  select(sname, ESA.Listing.Status, Entity.Description, ESA.Listing.Date, Common.Name) %>%
  as_tibble()

View(FWS_spp)

FWS_spp %>% count(ESA.Listing.Status) %>% arrange(desc(n))

FWS_spp %>% filter(ESA.Listing.Status!="Resolved Taxon") %>%
  mutate(snam_fix=recode(sname,
                         "Lithobates cascadae"= "Rana cascadae",
                         "Rana areolata sevosa"  = "Lithobates sevosus"  , 
                         "Rana yavapaiensis"= "Lithobates yavapaiensis",
                         "Rana subaquavocalis" = "Lithobates subaquavocalis",
                         "Rana tarahumarae" = "Lithobates tarahumarae",
                         "Rana aurora aurora" = "Rana aurora"
  )) %>% #fix some taxa names
  left_join(RCS_index, by=c('snam_fix'='species')) %>% 
  select(snam_fix, ESA.Listing.Status, RCS_WS) %>% View()

# look up these frogs here: https://ecos.fws.gov/ecp/report/species-candidate-removed-or-withdrawn
#also Rana subaquavocalis 

anuran_genus
NatGCC_list %>%
  filter(genus %in% anuran_genus) %>%
  arrange(States.2015)

# Figure 4: Genera description ----
ord_gen<-RCS_Data %>%
  mutate(genera=gsub("\\ .*",'',species)) %>%
  group_by(genera) %>% summarise(medRCS=median(RCS_WS)) %>%
  arrange(medRCS) %>% pull(genera)
taxa_table<-read.csv('gbif_anuran_taxa_info.csv')

library(scales)
oh<-NULL
gen_df<-RCS_Data %>% 
  mutate(genera=gsub("\\ .*",'',species)) %>%
  group_by(genera) %>%
  mutate(n_sp=n(),
         GenFac=factor(genera, levels=ord_gen)) %>%
  left_join(taxa_table) %>% 
  left_join(con_info, by=c('species'='scientific_name')) 
gen_con_df<-gen_df  %>%
  pivot_longer(cols=c("IUCNcategory", "ESA.Listing.Status",
                      "NatSGCN")) %>%
  count(name, n_sp, value)%>%
  group_by(genera, name, n_sp) %>%
  mutate(perc=n/n_sp) %>% 
  filter(value %in% c("Threatened","Endangered","VU","EN","in =1 SWAP")) %>%
  summarize(tot_p=paste0(round(sum(perc),2)*100, "%")) %>%
  pivot_wider(values_from=tot_p, values_fill="0%") %>%
  mutate(graph_text=paste0(n_sp, " (",IUCNcategory,"/", 
                           ESA.Listing.Status,"/",
                           NatSGCN, ")")) %>%
  bind_rows(data.frame(genera=c("Eleutherodactylus"),
                       n_sp=c(3),
                       graph_text=c("3 (0%/0%/0%)")))
con_bar<-gen_con_df %>%
  rename(genus=genera) %>% left_join(taxa_table %>% count(family, genus)) %>%
  pivot_longer(cols=c(NatSGCN,ESA.Listing.Status, IUCNcategory)) %>%
  select(genus, n_sp, name, value, family) %>%
  mutate(per=as.numeric(gsub('\\%','',value))/100,
         GenFac=factor(genus, levels=rev(ord_gen)),
         con.lab=factor(name, 
                        levels=c("IUCNcategory","ESA.Listing.Status","NatSGCN"),
                        labels=c("IUCN","ESA","SGCN")))%>%
  ggplot()+
  geom_col(aes(x=con.lab, y=per, fill=family))+
  facet_grid(GenFac~.)+
  scale_y_continuous('',labels=scales::percent, expand=c(0,0),
                     breaks=c(0,.5,1), position='right')+
  scale_x_discrete('')+
  scale_fill_viridis_d('Family')+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text=element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y.right = element_text(size=7),
        axis.text.x=element_text(size=8, angle = 15),
        legend.text = element_text(size=8),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-5,-5,-5,-5))
con_bar

gen_bars<-gen_df %>% left_join(con_plot_df, by=c('species'='scientific_name')) %>%
  ggplot()+
  geom_vline(xintercept=c(0.75,.5, .25), color="lightgrey")+
  geom_point(aes(y=GenFac, x=RCS_WS, color=`Listed on:`), 
             position=position_jitter(height=0.25),
             alpha=0.7)+
  geom_boxplot(aes(y=GenFac, x=RCS_WS), fill=NA,
               outlier.color = NA)+
  geom_rect(aes(fill=family, xmin=1.01, xmax=Inf, 
                ymin=as.numeric(GenFac)-.5, 
                ymax=as.numeric(GenFac)+.5))+
  geom_text(data=gen_con_df,
            aes(y=genera, x=.1, label=n_sp),
            size=3)+
  scale_x_reverse("Watershed RCS", limits=c(1.05,0), 
                  expand=expansion(mult = c(0.01, 0.01)),
                  breaks=c(1,.75,.5,.25,0))+
  scale_y_discrete("Genus")+
  scale_fill_viridis_d(guide=F)+
  scale_color_manual(values=c(viridis_pal(option='magma', 
                                          begin=.3, end=.75)(3), 
                              'lightgrey'),
                     guide=F)+
  #facet_wrap(~family, scales='free_y', ncol=1,
  #               strip.position='right')+
  theme_cowplot()+
  theme(axis.text.x = element_text(size=9),
        axis.text.y = element_text(size=9, face='italic'),
        axis.title.x = element_text(size=9),
        axis.title.y=element_text(size=9),
        strip.background = element_rect(fill=NA),
        strip.text=element_text(size=9, angle=90),
        strip.placement = 'outside',
        legend.position = 'bottom')

gen_bars
plot_grid(gen_bars, con_bar, rel_widths=c(1, .7), nrow=1)
ggsave("rcs_results/figures/Fig4_gen_wofam.jpeg", width=6.5, height=6)


#tutorial https://stackoverflow.com/questions/52341385/how-to-automatically-adjust-the-width-of-each-facet-for-facet-wrap
library(grid)
ohp = ggplot_gtable(ggplot_build(oh))
#gtable::gtable_show_layout(ohp)
# get the number of unique x-axis values per facet (1 & 3, in this case)
y.var <- sapply(ggplot_build(oh)$layout$panel_scales_y,
                function(l) length(l$range$range))
# change the relative widths of the facet columns based on
# how many unique x-axis values are in each facet
ohp$heights[ohp$layout$t[grepl("panel", ohp$layout$name)]] <- ohp$heights[ohp$layout$t[grepl("panel", ohp$layout$name)]] * y.var
grid.draw(ohp)
jpeg('rcs_results/figures/Fig3_genbox.jpg',
     width=3.5, height=7, units="in", res=150); grid.draw(ohp); dev.off()
#ggsave('rcs_results/figures/Fig3_genbox.jpg', height=6, width=3.3)

# Figure S4 area and state swaps -----
AOOs <- read.csv('rcs_results/intermediate results/AOO HUC12 Output_20210331.csv')

nswap_df<-con_info %>% 
  left_join(AOOs, by='scientific_name') %>%
  left_join(RCS_index_con, by=c("scientific_name"="species"))%>%
  filter(!is.na(RCS_WS),
         spatial_extent=='entire native NA range') %>%
  select(scientific_name, watershed, RCS_WS, nStates_2015)%>%
  replace_na(list('nStates_2015'=0)) %>%
  mutate(cl=ifelse(nStates_2015==0, "Not Listed", "Listed"))
library(ggrepel); library(scales)
ggplot(mapping=aes(x=nStates_2015, y=watershed))+
  geom_point(data=nswap_df, aes(size=RCS_WS, color=cl, shape=cl), alpha=0.5)+
  geom_text_repel(data=nswap_df %>% filter(nStates_2015 > 10),
                  aes(label=scientific_name), size=3)+
  scale_y_log10(expression("Area of Occupancy, km"^ 2),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  scale_x_continuous("Times listed on a SWAP")+
  scale_shape_manual(values=c(16,18), guide=F)+
  scale_color_manual(values=c("grey","black"), guide=F)+
  scale_size_continuous(guide=F)+
  theme_cowplot()
ggsave("rcs_results/figures/FigS4_swap_aoo.jpg", width=4, height=4)
nswap_df %>% filter(nStates_2015==0)

RCS_index %>% filter(!(species %in% nswap_df$final.taxa)) %>% pull(species) 
