library(raster)
setwd("W:/UR_intern/UR2/SW1-1/Mitarbeiter/Petermann/Paper/Mapping Indoor Rn - geogenic hazard and actual risk/GitHub/Data/")
#setwd("GitHub/")

#load predictor data
GRHI<-raster("GRHI.tif")
#load indoor Rn data
IRC<-read.csv("Indoor Rn.csv",sep=",")[,2:4]

#encode indoor Rn concentration exceedances
IRC$Exc_100<-as.factor(ifelse(IRC$Rn_cor>=100,1,0))
IRC$Exc_300<-as.factor(ifelse(IRC$Rn_cor>=300,1,0))

##########                          #########
###  Fit logistic regression models       ###  
#########                           #########

glm100.GRHI <- glm(Exc_100~GRHI,data=IRC,
                   family=binomial(link="logit"))

glm300.GRHI <- glm(Exc_300~GRHI,data=IRC,
                   family=binomial(link="logit"))


##########                          #########
###  Fit calibration model                ###  
#########                           #########

library(mgcv)
## 100 Bq/m³
#calculate probabilities of logistic regression model
fitted.response100<-predict(glm100.GRHI,type="response")
#fit ordinary GAM for smoothing
cal.mod.GRHI_100<-gam(IRC$Exc_100~s(fitted.response100,bs="ts",m=3),family="binomial") #s() defines smooths in GAM formula; bs="ts" thin plate regression splines defines smoothing basis

## 300 Bq/m³
fitted.response300<-predict(glm300.GRHI,type="response")
#fit ordinary GAM for smoothing
cal.mod.GRHI_300<-gam(IRC$Exc_300~s(fitted.response300,bs="ts",m=3),family="binomial") #s() defines smooths in GAM formula; bs="ts" thin plate regression splines defines smoothing basis


#######   ######
###  Mapping ###  
######    ######

#define grid
x=seq(280500,921500,1000)
y=seq(5235500,6106500,1000)
xy=expand.grid(x,y)
#coordinate reference system UTM 32
UTM32=CRS("+init=epsg:25832")


### logistic model

# generate estimation grid
Grid=SpatialPointsDataFrame(xy[,1:2],xy,proj4string = UTM32,bbox=NULL)
# extract GRHI data for estimation grid
Grid@data$GRHI<-extract(GRHI,Grid)
# name variable as predictor
Grid<-as.data.frame(Grid["GRHI"])

#####100 Bq/m³

#make glm prediction
Prob.raw100 <- predict(glm100.GRHI,newdata=Grid["GRHI"],type="response")

# calibration
rawdata100<-as.data.frame(Prob.raw100)
colnames(rawdata100)="fitted.response100"
Prob100.cal<-predict(cal.mod.GRHI_100,newdata=rawdata100,type="response")

# generate sp object
Prob100.cal.p<-SpatialPointsDataFrame(xy[,1:2],as.data.frame(Prob100.cal),proj4string = UTM32,bbox=NULL)
# grid sp object
gridded(Prob100.cal.p) <- TRUE
# generate raster
Prob100.cal.r <- raster(Prob100.cal.p)


#####300 Bq/m³
### logistic model

#make glm prediction
Prob.raw300 <- predict(glm300.GRHI,newdata=Grid["GRHI"],type="response")

# calibration
rawdata300<-as.data.frame(Prob.raw300)
colnames(rawdata300)="fitted.response300"
Prob300.cal<-predict(cal.mod.GRHI_300,newdata=rawdata300,type="response")

# generate sp object
Prob300.cal.p<-SpatialPointsDataFrame(xy[,1:2],as.data.frame(Prob300.cal),proj4string = UTM32,bbox=NULL)
# grid sp object
gridded(Prob300.cal.p) <- TRUE
# generate raster
Prob300.cal.r <- raster(Prob300.cal.p)

setwd("W:/UR_intern/UR2/SW1-1/Mitarbeiter/Petermann/Paper/Mapping Indoor Rn - geogenic hazard and actual risk/GitHub/Data/results/")
# writeRaster(Prob100.cal.r,"Prob100_cal.tif")
# writeRaster(Prob300.cal.r,"Prob300_cal.tif")


#######   ######
###  Fig 6   ###  
######    ######

# install.packages("tmap")
# install.packages("tmaptools")
library(tmap)
library(tmaptools)

# define values breaks
breaks100 = c(0,.10,.20,.30,.40,.50,.60,Inf)
breaks300 = c(0,0.01,.025,.05,.075,.1,.15,Inf)

setwd("W:/UR_intern/UR2/SW1-1/Mitarbeiter/Petermann/Paper/Mapping Indoor Rn - geogenic hazard and actual risk/GitHub/Data/federal states administrative boundary/")
BL <- read_shape("VG250_LAN.shp",as.sf=TRUE) # bundeslaender -> federal states

#create tmap_plots
P100_GLM<-  tm_shape(BL)+
  tm_borders(col="black",lwd=1,alpha=0.8)+
  tm_shape(Prob100.cal.r*100,bbox=extent(BL))+
  tm_raster(palette = "Reds", breaks = breaks100*100,style="fixed",title="[%]")+
  tm_layout(legend.show=TRUE,
            panel.labels = "Probability IRC>100 Bq/m³",panel.label.size = 1,legend.outside=FALSE,
            legend.position=c("LEFT","TOP"))+
  tm_shape(BL)+
  tm_borders(col="black",lwd=1,alpha=0.8)+
  tm_scale_bar(position = c("RIGHT", "BOTTOM"),width=0.15)+
  tm_credits("(c) GeoBasis-DE / BKG 2020", position=c("LEFT", "BOTTOM"))

P300_GLM <- tm_shape(BL)+
  tm_borders(col="black",lwd=1,alpha=0.8)+
  tm_shape(Prob300.cal.r*100,bbox=extent(BL))+
  tm_raster(palette = "Reds", breaks = breaks300*100,style="fixed",title="[%]")+
  tm_layout(legend.show=TRUE,
            panel.labels = "Probability IRC>300 Bq/m³",panel.label.size = 1,legend.outside=FALSE,
            legend.position=c("LEFT","TOP"))+
  tm_shape(BL)+
  tm_borders(col="black",lwd=1,alpha=0.8)+
  tm_scale_bar(position = c("RIGHT", "BOTTOM"),width=0.15)+
  tm_credits("(c) GeoBasis-DE / BKG 2020", position=c("LEFT", "BOTTOM"))

#plot
tmap_arrange(P100_GLM,P300_GLM,ncol=2,nrow=1)

#######   ######
###  Fig 9   ###  
######    ######

breaksRPA10=c(0,.1,1)
breaksRPA5=c(0,.05,1)
P300_GLM_RPA10=  tm_shape(BL)+
  tm_borders(col="black",lwd=1,alpha=0.8)+
  tm_shape(Prob300.cal.r)+
  tm_raster(palette = "Reds", breaks = breaksRPA10,style="fixed",title="Probability")+
  tm_layout(legend.show=TRUE,
            panel.labels = "Probability(IRC>300 Bq/m³) > 10 %",panel.label.size = 1,legend.outside=FALSE,
            legend.outside.position="bottom",legend.position=c("left","top"))+
  tm_shape(BL)+
  tm_borders(col="black",lwd=1,alpha=0.8)+
  tm_scale_bar(position = c("RIGHT", "BOTTOM"),width=0.15)+
  tm_credits("(c) GeoBasis-DE / BKG 2020", position=c("LEFT", "BOTTOM"))
P300_GLM_RPA5=  tm_shape(BL)+
  tm_borders(col="black",lwd=1,alpha=0.8)+
  tm_shape(Prob300.cal.r)+
  tm_raster(palette = "Reds", breaks = breaksRPA5,style="fixed",title="Probability")+
  tm_layout(legend.show=TRUE,
            panel.labels = "Probability(IRC>300 Bq/m³) > 5 %",panel.label.size = 1,legend.outside=FALSE,
            legend.outside.position="bottom",legend.position=c("left","top"))+
  tm_shape(BL)+
  tm_borders(col="black",lwd=1,alpha=0.8)+
  tm_scale_bar(position = c("RIGHT", "BOTTOM"),width=0.15)+
  tm_credits("(c) GeoBasis-DE / BKG 2020", position=c("LEFT", "BOTTOM"))

setwd("W:/UR_intern/UR2/SW1-1/Mitarbeiter/Petermann/Projekte/Hazard Index - Risk - Exceedence Prob/Figures/")
#tiff("RPA 10 vs. 5.tif",width=22,height=13,units="cm",res=150)
tmap_arrange(P300_GLM_RPA10,P300_GLM_RPA5, ncol=2,nrow=1)
#dev.off()


rawdata300<-as.data.frame(c(0.002))
colnames(rawdata300)="fitted.response300"
Prob300.cal<-predict(cal.mod.GRHI_300,newdata=rawdata300,type="response")
