library(tidyverse) 
amphibio<-read.csv('C:/Users/Owner/Downloads/AmphiBIO_v1/AmphiBIO_v1.csv')
names(amphibio)
head(amphibio)
amphibio %>% as_tibble() %>% 
  filter(Order=="Anura") %>%
  group_by(Order, Family, Species) %>%
  tally()
tax_reference<-read.csv("data/AnuranTaxRef_20200708.csv")
tax_reference
amphibio_anurans<-amphibio %>% as_tibble() %>% 
  filter(Order=="Anura")
# checking the taxa I use are all valid
library(ritis)
testTax<-function(taxa){
  vec<-NULL
  result<-NULL
  for(w in 1:length(taxa)){
    result<-itis_search(q=paste0('nameWOInd:',gsub(' ','\\ ', taxa[w], fixed=T)))$usage
    if(is.null(result)){
      vec[w]<-1}else{
      vec[w]<-length(grep('invalid',result))}
  }
  print(paste('A total of ', length(which(vec==1)), 'of', length(taxa), 'are invalid or missing.'))
  return(amphibio_anurans$Species[which(vec==1)])
}

res<-testTax(amphibio_anurans$Species)

#does res match any of other invalid taxa?
filter(tax_reference, tax_reference$species %in% res)


amphibio_anurans_filt<- amphibio_anurans %>%
  filter(Species %in% tax_reference$species) %>%
  left_join(tax_reference, by=c("Species"="species"))
amphibio_anurans_filt %>% 
  group_by(final.taxa) %>% 
  summarize(n_obs=n()) %>%
  arrange(desc(n_obs))
nrow(amphibio_anurans_filt)
length(unique(amphibio_anurans_filt$final.taxa))
# so 103 taxa in this database
# need to investigate how they built taxa table

View(amphibio_anurans_filt)
library(tidyverse)
globalTherm<-read.csv('C:/Users/Owner/Downloads/GlobTherm/GlobalTherm_upload_02_11_17.csv')
names(globalTherm)
anuran_therm<-globalTherm %>% as_tibble() %>%
  filter(Order=="Anura") %>%
  mutate(taxa=paste(Genus, Species)) %>%
  dplyr::select(taxa, everything())
anuran_therm %>%
  #filter(taxa %in% tax_reference$final.taxa) %>%
  group_by(Genus, max_metric) %>% tally()
anuran_therm %>% 
  filter(taxa %in% tax_reference$final.taxa) %>%
  dplyr::select(taxa, max_metric, lat_min, long_min, location_min, elevation_min)


#### what is in the literature? ####
library(rwos)
sid <- wos_authenticate() #have to be at a campus (or VPN) that has access to WoS
res <- wos_search(sid, "TS=Anura* tadpole temperature min*")
pubs <- wos_retrieve(res)
View(pubs)
