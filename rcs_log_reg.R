# RCS to predict whether considered threatened or not ----
# Last edited by Traci Dubose ~ 2020Oct15
#stat notes & examples
# using Practical Guide to Logistic Regression by Joseph M. Hilbe (2015) as a guide
# also https://m-clark.github.io/docs/GAMS.pdf for explanation of GAMS
library(mgcv)

# Split data and make conservation status a binary variable ----
rcs_na<-RCS_index_con %>%
  mutate(ESA_bin=case_when(`Endangered Species Act`=='Not on ESA'~0,
                           T~1),
         State_bin=case_when(`Species of Greatest\nConserv. Need`=='in >1 SWAP'~1,
                             T~0),
         Int_bin=case_when(`IUCN Red List` %in% c('Not listed', 'Least Concern', 'Near Threatened')~0,
                           T~1)) %>%
  filter(spatial_extent=='entire native NA range')
rcs_l48<-RCS_index_con %>%
  mutate(ESA_bin=case_when(`Endangered Species Act`=='Not on ESA'~0,
                           T~1),
         State_bin=case_when(`Species of Greatest\nConserv. Need`=='in >1 SWAP'~1,
                             T~0),
         Int_bin=case_when(`IUCN Red List` %in% c('Not listed', 'Least Concern', 'Near Threatened')~0,
                           T~1)) %>%
  filter(spatial_extent=='native L48 range')
hist(rcs$RCS_WS)

# Run logistic regression models ----
esagam_er <- gam(ESA_bin~s(RCS_WS), family=binomial, data=rcs_na)
esaglm_er <- glm(ESA_bin~RCS_WS, family=binomial, data=rcs_na)
summary(esaglm)
AIC(esagam_er, esaglm_er)
esagam_l4 <- gam(ESA_bin~s(RCS_WS), family=binomial, data=rcs_l48)
esaglm_l4 <- glm(ESA_bin~RCS_WS, family=binomial, data=rcs_l48)
AIC(esagam_l4, esaglm_l4)
deviance(esagam_er); deviance(esagam_l4)

intgam_er <- gam(Int_bin~s(RCS_WS), family=binomial, data=rcs_na)
intgam_l4 <- gam(Int_bin~s(RCS_WS), family=binomial, data=rcs_l48)
deviance(intgam_er); deviance(intgam_l4)

stagam_er <- gam(State_bin~s(RCS_WS), family=binomial, data=rcs_na)
stagam_l4 <- gam(State_bin~s(RCS_WS), family=binomial, data=rcs_l48)
deviance(stagam_er); deviance(stagam_l4)

# Outcome table
log.p.vals<-rbind(summary(intgam_er)$s.table,
                  summary(esagam_er)$s.table,
                  summary(stagam_er)$s.table) %>%
  as_tibble() %>%
  mutate(name=c('IUCN','ESA','SWAP'))


# Plot the model deviance for goodness of fit ----
log.reg.dev.p<-data.frame(model=c("esagam_er","esagam_l4","intgam_er","intgam_l4","stagam_er","stagam_l4"),
           deviance=c(deviance(esagam_er), deviance(esagam_l4), deviance(intgam_er), deviance(intgam_l4),deviance(stagam_er), deviance(stagam_l4)),
           con_org=factor(c("ESA","ESA","IUCN","IUCN","SGCN","SGCN"),
                          levels=c("SGCN","ESA","IUCN")),
           spat_extent=factor(rep(c("entire","l48"), 3),
                              levels=c("l48","entire"))) %>%
  ggplot()+
  geom_col(aes(x=spat_extent, y=deviance, fill=model))+
  facet_wrap(~con_org, scales="free_x")+
  ggtitle('Model Fit')+
  scale_fill_viridis_d('List', option='E', guide=F)+
  scale_y_continuous('Model Deviance', limits=c(0,100), expand = c(0,0))+
  scale_x_discrete('Model spatial extent')+
  theme_cowplot()
log.reg.predi.p<-data.frame(model=c("esagam_er","esagam_l4","intgam_er","intgam_l4","stagam_er","stagam_l4"),
           deviance=c(60,65,50,55,70,75),
           con_org=factor(c("ESA","ESA","IUCN","IUCN","SGCN","SGCN"),
                          levels=c("SGCN","ESA","IUCN")),
           spat_extent=factor(rep(c("entire","l48"), 3),
                              levels=c("l48","entire"))) %>%
  ggplot()+
  geom_col(aes(x=spat_extent, y=deviance, fill=model))+
  facet_wrap(~con_org, scales="free_x")+
  ggtitle('Predictions')+
  scale_y_continuous('Model Deviance')+
  scale_x_discrete('Model spatial extent')+
  theme_cowplot()+
  scale_fill_viridis_d('List', option='E', guide=F)+
  theme(axis.text.y = element_text(color=NA))

# Plot logistic regression predictions ----
new_data_fake<-data.frame(State_bin=1,
                          ESA_bin=1,
                          Int_bin=1,
                          RCS_WS=seq(0,1,by=0.01))
predicted_data<-new_data_fake %>%
  mutate('IUCN_er'=predict(intgam_er, newdata=new_data_fake,type="response"),
         'ESA_er'=predict(esagam_er, newdata=new_data_fake,type="response"),
         'SWAP_er'=predict(stagam_er, newdata=new_data_fake,type="response"),
         'IUCN_l4'=predict(intgam_l4, newdata=new_data_fake,type="response"),
         'ESA_l4'=predict(esagam_l4, newdata=new_data_fake,type="response"),
         'SWAP_l4'=predict(stagam_l4, newdata=new_data_fake,type="response"))%>%
  pivot_longer(cols = -c('State_bin',"ESA_bin", "Int_bin","RCS_WS"))

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
  rename('ESA'=ESA_bin, 'SWAP'=State_bin, 'IUCN'=Int_bin) %>%
  pivot_longer(cols = c('ESA',"SWAP", "IUCN")) %>% #View()%>%
  ggplot()+
  geom_line(data=predicted_data %>% rename(model='name') %>%
              mutate(name=sub("_.*","",model)), 
            aes(x=RCS_WS, y=value, 
                color=model), size=1.5)+
  geom_point(aes(x=RCS_WS, y=value, 
                                 fill=as.factor(value)),
             alpha=0.3, position=position_jitter(height=.02),
             shape=21)+
  #geom_text(data=log.p.vals, x=0, y=1,
  #          aes(label=paste('p = ', round(`p-value`,2))),
  #          hjust=0)+
  scale_x_continuous('watershed RCS')+
  scale_y_continuous('Prob. of being listed')+
  
  #scale_linetype_manual('List',values=c(1,1,2,2,4,4))+
  scale_colour_viridis_d('List', option='E', guide=F)+
  scale_fill_manual('Listed as\nThreatened',values=c('white', 'black'),
                    labels=c('No','Yes'))+
  facet_wrap(~name, ncol=1)+
  theme_cowplot()+  theme(legend.position = 'bottom',
                          axis.text.x = element_text(size=10),
                          axis.text.y = element_text(size=10),
                          legend.justification = 'center') +
  guides(fill=guide_legend(nrow=2,
                           override.aes=list(alpha=1)))
ggsave('rcs_results/figures/logreg_output.jpg', width=3.5, height=7)

# Figure XXX ----
plot_grid(tall.log.pred.p,
          plot_grid(log.reg.predi.p, log.reg.dev.p, ncol=1),
          rel_widths = c(1,.75))
ggsave('rcs_results/figures/logreg_pred.jpg', width=6.5, height=7)


# Old code ------
summary(esagam_er)
summary(intgam)
summary(stagam)

str(intgam)
log_res_df<-rcs %>% as_tibble() %>%
  mutate(IUCNres=intgam$residuals,
         ESAres=esagam$residuals,
         Stateres=stagam$residuals)
  

library(ggrepel)
ip<-ggplot(data=log_res_df,
           aes(x=spatial_extent, y=IUCNres))+
  geom_point(alpha=0.2)+
  geom_text_repel(data=log_res_df[log_res_df$IUCNres>3,],
                  aes(label=species),
                  size=3)+
  theme_cowplot()



ep<-ggplot(data=log_res_df)+
  stat_summary(aes(x=ESA_bin, y=ESAres))+
  scale_y_continuous(lim=c(-4, 9))+
  facet_wrap(~spatial_extent)
spp<-ggplot(data=log_res_df)+
  stat_summary(aes(x=State_bin, y=Stateres))+
  scale_y_continuous(lim=c(-4, 9))+
  facet_wrap(~spatial_extent)

plot_grid(ip,ep, spp, ncol=1)

set.seed(2) ## simulate some data... 
dat <- gamSim(1,n=400,dist="normal",scale=2)
b <- gam(y~s(x0)+s(x1)+s(x2)+s(x3),data=dat)
summary(b)

intgam
