library(rgdal)
library(tmap)
library(tmaptools)
library(sf)
setwd("W:/UR_intern/UR2/SW1-1/Mitarbeiter/Petermann/Projekte/Hazard Index - Risk - Exceedence Prob/GIS/")
buildings<-readOGR(".","Municipalities with buildings exc100 300 final")
buildings@data$p300_mean<-buildings@data$p300_mean*100
buildings@data$p100_mean<-buildings@data$p100_mean*100
breaks_build100=c(0,50,200,500,2000,5000,Inf)
breaks_build300=c(0,10,50,200,500,2000,Inf)
breaks100 = c(0,.10,.20,.30,.40,.50,.60,Inf)
breaks300 = c(0,0.01,.025,.05,.075,.1,.15,Inf)

setwd("W:/UR_intern/UR2/SW1-1/Mitarbeiter/Petermann/Daten/GIS/")
UTM32=CRS("+init=epsg:25832")
BL <- read_shape("VG250_LAN.shp",as.sf=TRUE)

probs100=  tm_shape(BL)+
  tm_fill(col="lightgray")+
  tm_borders(col="black",lwd=1,alpha=0.8)+
  tm_shape(buildings)+
  tm_fill(palette = "Reds",col="p100_mean", breaks = breaks100*100,style="fixed",title="[%]")+
  tm_layout(legend.show=TRUE,
            panel.labels = "Probability (IRC>100 Bq/m³)",panel.label.size = 1,legend.outside=FALSE,
            legend.outside.position="bottom",legend.position=c("LEFT","TOP"))+  
  tm_shape(BL)+
  tm_borders(col="black",lwd=1,alpha=0.8)+
  tm_scale_bar(position = c("RIGHT", "BOTTOM"),width=0.15)+
  tm_credits("(c) GeoBasis-DE / BKG 2020", position=c("LEFT", "BOTTOM"))


probs300=  tm_shape(BL)+
  tm_fill(col="lightgray")+
  tm_borders(col="black",lwd=1,alpha=0.8)+
  tm_shape(buildings)+
  tm_fill(palette = "Reds",col="p300_mean", breaks = breaks300*100,style="fixed",title="[%]")+
  tm_layout(legend.show=TRUE,
            panel.labels = "Probability (IRC>300 Bq/m³)",panel.label.size = 1,legend.outside=FALSE,
            legend.outside.position="bottom",legend.position=c("LEFT","TOP"))+
  tm_shape(BL)+
  tm_borders(col="black",lwd=1,alpha=0.8)+
  tm_scale_bar(position = c("RIGHT", "BOTTOM"),width=0.15)+
  tm_credits("(c) GeoBasis-DE / BKG 2020", position=c("LEFT", "BOTTOM"))


buildings100=  tm_shape(BL)+
  tm_fill(col="lightgray")+
  tm_borders(col="black",lwd=1,alpha=0.8)+
  tm_shape(buildings)+
  tm_fill(palette = "Reds",col="geb_p100", breaks = breaks_build100,style="fixed",title="Buildings")+
  tm_layout(legend.show=TRUE,
            panel.labels = "Residential buildings >100 Bq/m³",panel.label.size = 1,legend.outside=FALSE,
            legend.outside.position="bottom",legend.position=c("LEFT","TOP"))+
  tm_shape(BL)+
  tm_borders(col="black",lwd=1,alpha=0.8)+
  tm_scale_bar(position = c("RIGHT", "BOTTOM"),width=0.15)+
  tm_credits("(c) GeoBasis-DE / BKG 2020", position=c("LEFT", "BOTTOM"))


buildings300=  tm_shape(BL)+
  tm_fill(col="lightgray")+
  tm_borders(col="black",lwd=1,alpha=0.8)+
  tm_shape(buildings)+
  tm_fill(palette = "Reds",col="geb_p300", breaks = breaks_build300,style="fixed",title="Buildings")+
  tm_layout(legend.show=TRUE,
            panel.labels = "Residential buildings >300 Bq/m³",panel.label.size = 1,legend.outside=FALSE,
            legend.outside.position="bottom",legend.position=c("LEFT","TOP"))+
  tm_shape(BL)+
  tm_borders(col="black",lwd=1,alpha=0.8)+
  tm_scale_bar(position = c("RIGHT", "BOTTOM"),width=0.15)+
  tm_credits("(c) GeoBasis-DE / BKG 2020", position=c("LEFT", "BOTTOM"))


setwd("W:/UR_intern/UR2/SW1-1/Mitarbeiter/Petermann/Projekte/Hazard Index - Risk - Exceedence Prob/Figures/")
tiff("Buidlings above 100_300 Bqm3_.tif",width=26,height=26,units="cm",res=150)
tmap_arrange(probs100,probs300,buildings100,buildings300, ncol=2,nrow=2)
dev.off()
