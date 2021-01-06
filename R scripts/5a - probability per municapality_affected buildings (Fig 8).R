install.packages("spatialEco")
install.packages("rgdal")
library(raster)
library(tmap)
library(tmaptools)
library(rgdal)

#coordinate reference system UTM 32
UTM32=CRS("+init=epsg:25832")

setwd("W:/UR_intern/UR2/SW1-1/Mitarbeiter/Petermann/Paper/Mapping Indoor Rn - geogenic hazard and actual risk/GitHub/Data/results/")
Prob100_cal<-raster("Prob100_cal.tif",crs=UTM32)
Prob300_cal<-raster("Prob300_cal.tif")

setwd("W:/UR_intern/UR2/SW1-1/Mitarbeiter/Petermann/Paper/Mapping Indoor Rn - geogenic hazard and actual risk/GitHub/Data/built-up area/")
builtup_area <- readOGR(".","built-up area") # bundeslaender -> federal states


builtup_area<-spTransform(builtup_area,UTM32)

#dev.off()
library(rgdal)
setwd("W:/UR_intern/UR2/SW1-1/Mitarbeiter/Petermann/Projekte/Hazard Index - Risk - Exceedence Prob/GIS/")
builtup<-readOGR(".","settlements areas incl buildings")
setwd("W:/UR_intern/UR2/SW1-1/Mitarbeiter/Petermann/Daten/Bundesamt für Kartographie und Geodäsie/Verwaltungsgrenzen Deutschland/vg250_ebenen/")
munic<-readOGR(".","VG250_GEM")

library(dplyr)
joined_tibble <- right_join(builtup@data,munic@data,by="AGS")
municip_exc <- merge(x=munic,builtup@data,by="AGS",all.x=TRUE)
municip_exc@data$geb_p100 <- as.numeric(as.character(municip_exc@data$geb_p100))
municip_exc@data$geb_p300 <- as.numeric(as.character(municip_exc@data$geb_p300))
setwd("W:/UR_intern/UR2/SW1-1/Mitarbeiter/Petermann/Projekte/Hazard Index - Risk - Exceedence Prob/GIS/")
#writeOGR(municip_exc[,c("AGS","GEN.x","Wohngbd","p300_mean","p100_mean","geb_p100","geb_p300")],dsn=".","Municipalities with buildings exc100 300",driver="ESRI Shapefile",overwrite=TRUE)
