#TPD early notes for anuran rcs product for this summer

#looking at some of the data
library(tidyverse); library(sf); library(rgeos); library(sp); library(rgdal)
library(cowplot)

#ideas - chunk years - start at ~30 year increments? 30 - 20 - 10 
# i don't think 5 is good. its about a funding cycle (in my brain).

#sample size affect of rarity
# .25, .5, .7, .8, .9, .1? 
# or threshold? 

#how is sample size correlated with time? What time periods reduce this correlation?

list.files('G:/Shared drives/Anuran Biodiversity Project/SDMs/Occurrence Data/Download & Filtering/')[1:4]

for(k in c(1:4)){
  #as each sp has its own folder, this pulls out sp name at position k
  a_sp<-list.files('G:/Shared drives/Anuran Biodiversity Project/SDMs/Occurrence Data/Download & Filtering/')[k]
  #writing the correct file name to pull the csv file out of
  file_a_sp<-paste('G:/Shared drives/Anuran Biodiversity Project/SDMs/Occurrence Data/Download & Filtering/',
               a_sp, sep='')
  #assigning each csv to the correct G.species object
  assign(paste(substr(a_sp,1,1), #pull genus inital out
               gsub("[A-z ]* ","", a_sp), #keep everything after the first space
               sep='.'),
         read.csv(paste(file_a_sp, '/', list.files(file_a_sp, pattern = 'v2.csv'), sep='')))#read the file that has v2.csv in the file name
  #check this actually worked correctly
  print(paste(paste(substr(a_sp,1,1), gsub("[A-z ]* ","", a_sp), sep='.'), #object name
              unique(get(paste(substr(a_sp,1,1), gsub("[A-z ]* ","", a_sp),#species within that object name
                               sep='.'))$species), sep="=="))

  }



#sp_5 doesn't have a merged file yet?

graphdf<-bind_rows(mget(paste('sp',c(1:4,6:10), sep='_')))

min(graphdf$year);max(graphdf$year)
makeBreaks<-function(df, period){
  df_all<-NULL
  for(i in unique(df$species)){
    df_sp <- df[df$species==i,]
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
year_counts<-makeBreaks(graphdf, c(10,20,30,50))

ggplot(year_counts) +
  geom_point(aes(x=mids, y=n, color=as.factor(period)), alpha=0.5)+
  facet_wrap(~species, scales = "free")+
  scale_color_discrete(name="period")+
  theme_bw()


one_sp<-graphdf[graphdf$species=='Acris blanchardi',]
hist(one_sp$year, breaks)


applyBreaks<-function(df, period){
  perds<-makeBreaks(df, period)[,c("period","min","mids","max")] %>%
    filter(!duplicated(.))
  df[,"mids"]<-perds$mids[findInterval(df$year, perds$min)]
  new_df<-df %>% left_join(perds)
  return(new_df)
}

subsample<-graphdf %>% sample_n(100)
t<-applyBreaks(subsample, c(50))  
which(!(t$year>=t$min));which(!(t$year<=t$max)) #want two integer(0)

#next step apply RCS at these time periods
# calculate rarity & climate sensitivity using the time periods within min-max

rarity<-function(df, period){
  new_df<-applyBreaks(df, period)
  for(q in period){
    p_df<-new_df %>% filter(period==q)
    p_df_sf<-st_as_sf(p_df,
                      coords=c("Longitude","Latitude"),
                      crs="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs")
  }
  return(p_df_sf)
}

plot(rarity(one_sp, 50)["year"])

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


