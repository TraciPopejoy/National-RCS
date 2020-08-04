# scripts/notes to investigate the difference between prism and worldclim data
library(tidyverse)
library(raster)

?getData

tax_reference<-read.csv("data/AnuranTaxRef_20200708.csv")
anuran.taxa<-unique(tax_reference$final.taxa)
head(anuran.taxa)

geodata<-read.csv("data/occ_data_used/Acris blanchardi.csv") %>%
  dplyr::select(family, genus, species, Longitude, Latitude, year, source, final.taxa)
geodata.sp<-SpatialPointsDataFrame(coords= cbind(geodata$Longitude, geodata$Latitude),
                                   data= geodata,
                                   proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))
head(geodata)

test<-getData('worldclim', var='bio', res=0.5, 
        lon=geodata$Longitude[1], lat=geodata$Latitude[1])
test
raster::extract(test, geodata.sp[1,])
