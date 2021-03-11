# Thoughts on rcs through time
# late revised 2020July01 by tpd

library(tidyverse); library(lubridate)

Anurans_df<-read.csv("data/anuran_occ_all20200708.csv")
tax_reference<-read.csv("data/AnuranTaxRef_20200708.csv")

names(Anurans_df)
Anurans_df %>% group_by(source) %>% tally()
Anurans_df %>% filter(is.na(Date), source=="HM") %>% select(final.taxa, Date, source,year)
Anurans_df %>% filter(is.na(year)) %>% group_by(source) %>% tally()

Anurans_df_year<-Anurans_df %>% 
  filter(!is.na(year),  
         year < year('2020-06-01')) #I needed to add a year to files that had no year
#data originally imported in december 2019, so anything after is from me!

min(Anurans_df_year$year);max(Anurans_df_year$year)

makeBreaks<-function(df, period){
  df_all<-NULL
  for(i in unique(df$final.taxa)){
    df_sp <- df[df$final.taxa==i,]
    for(k in period){
      j<-seq(to=2020, from=1700, by=k)
      if(sum(j==2020)<1) j<-c(j,2020)
      h<-hist(df_sp$year, breaks = j, plot=F)
      dfx<-data.frame(species=i,
                      period=k,
                      min=h$breaks[-length(h$breaks)],
                      mids=h$mid,
                      max=h$breaks[-1],
                      n=h$counts,
                      dens=h$density)
      df_all<-bind_rows(df_all, dfx)
    }
  }
  return(df_all)
}
year_counts<-makeBreaks(Anurans_df_year, c(50))

pdf(paste0('rcs_results/yearly_counts_sp',
           format(Sys.Date(), '%Y%m%d'),'.pdf'))
ggplot(year_counts) +
  geom_point(aes(x=mids, y=n, color=as.factor(period)), alpha=0.5)+
  facet_wrap(~species, scales = "free")+
  scale_color_discrete(name="period")+
  theme_bw()
dev.off()

year_counts_plot<-year_counts %>% 
  mutate(genus=gsub(" .*$","",species)) %>%
  filter(min>1850)
length(unique(year_counts$genus))

year_counts %>% group_by(min,mids, max) %>%
  tally()

ggplot()+
  geom_jitter(data=year_counts[year_counts$n!=0,], 
              aes(x=mids, y=n), alpha=0.2,
              width=10)+
  theme_bw()+
  scale_y_log10("Number of Observations")+
  scale_x_continuous(name="Year of Obs.\nmidpoints of 50y periods",
    breaks=c(1725,1775,1825,1875, 1925, 1975, 2010))
ggsave('rcs_results/figures/sample_size_time_50y.jpg',
       width=4, height=4)

ggplot()+
  geom_density(data=Anurans_df_year, aes(x=year))+
  theme_bw()+
  scale_x_continuous(name= "Year of Observation")
ggsave('rcs_results/figures/sample_size_time.jpg', 
       width=6, height=4)

applyBreaks<-function(df, period){
     perds<-makeBreaks(df, period)[,c("period","min","mids","max")] %>%
       filter(!duplicated(.))
     df[,"mids"]<-perds$mids[findInterval(df$year, perds$min)]
     new_df<-df %>% left_join(perds, by="mids", na_matches="never")
     return(new_df)
}

testDF<-applyBreaks(Anurans_df_year, c(50))  
which(!(testDF$year>=testDF$min))
which(!(testDF$year<=testDF$max)) #want two integer(0)

#next step apply RCS at these time periods
# calculate rarity & climate sensitivity using the time periods within min-max


#### Seeing if counts increase through time ####
#tutorial; https://stats.idre.ucla.edu/r/dae/poisson-regression/
library(sandwich)
MultSp<-makeBreaks(graphdf, 30)
m1<-glm(n~mids:species, data=MultSp, family="poisson")
summary(m1)

cov.m1<-vcovHC(m1, type="HC0") #Cameron and Trivedi (2009)
std.err<-sqrt(diag(cov.m1))

slopedf <- as_tibble(cbind(Estimate= coef(m1), "Robust SE" = std.err,
                           "Pr(>|z|)" = 2 * pnorm(abs(coef(m1)/std.err), lower.tail=FALSE),
                           LL = coef(m1) - 1.96 * std.err,
                           UL = coef(m1) + 1.96 * std.err)) %>%
  mutate(species= sub('mids:species','', names(coef(m1)))) %>%
  filter(species!="(Intercept)")

ggplot(slopedf) + 
  geom_hline(yintercept=0, size=2, color='red')+
  geom_linerange(aes(x=species, ymax=UL, ymin=LL))+
  geom_point(aes(x=species, y=Estimate), size=2)+
  scale_y_continuous(name="slope estimate")+
  theme_cowplot()+ 
  theme(axis.text.x = element_text(angle=30, hjust=1))


