# install.packages("raster")
# install.packages("Smisc")
library(raster)
library(Smisc)

##########      #########
###  load data  ###  
#########       #########

setwd("W:/UR_intern/UR2/SW1-1/Mitarbeiter/Petermann/Paper/Mapping Indoor Rn - geogenic hazard and actual risk/GitHub/Data/")
#setwd("GitHub/")

#load predictor data
GRP<-raster("GRP.tif")
GRHI<-raster("GRHI.tif")
#load indoor Rn data
IRC<-read.csv("Indoor Rn.csv",sep=",")[,2:4]

#encode indoor Rn concentration exceedances
IRC$Exc_100<-as.factor(ifelse(IRC$Rn_cor>=100,1,0))
IRC$Exc_300<-as.factor(ifelse(IRC$Rn_cor>=300,1,0))

##########                                     #########
###  spatial cross-validated performance evaluation  ###  
#########                                      #########


### spatial block creation (coordinates of observations required --> sf or SpatialPoints object)

# install.packages("blockCV")
# library(blockCV)
# spatial_blocks <-list()
# for (j in 1:20){                            # j determines number of repitions
#   set.seed(j)
#   sb <- spatialBlock(speciesData = IRC,     # data
#                      species="Exc_300",     # variable of interest
#                      rasterLayer = GRP,     # background map
#                      theRange = 40000,      # size of the blocks in metres
#                      k = 5,                 # number of folds
#                      selection = "random",
#                      iteration = 100,       # find evenly dispersed folds
#                      biomod2Format = TRUE,
#                      xOffset = 0,           # shift the blocks horizontally: no
#                      yOffset = 0)           # shift the blocks vertically: no
#   for (k in 1:5){                           # k=5->for each fold
#     ind<-sb[[3]][,k]                        # [[3]]access biomodtable; [,k] access fold k
#     ind_train<-which(ind == TRUE)           # extract indices for training data
#     spatial_blocks[[(j-1)*5+k]]<-ind_train  # create a list with j*k spatial blocks with indices of training folds
#   }
# }

#loading spatial block data: 20 times repeated 5-fold spatially blocked data->100 spatial blocks with 40 km block size
sb100<-loadObject("spatialBlocks IRC Exc100.R")
sb300<-loadObject("spatialBlocks IRC Exc300.R")

# define ranges for bins
bin.thresh300 <- c(0,0.02,0.04,0.06,0.08,0.1,0.15,1) # 300 Bq/m³ exceedance
bin.thresh100 <- c(0,0.1,0.15,0.2,0.25,0.3,0.45,1)   # 100 Bq/m³ exceedance

# create matrices for storing performance results

# defining metrics to be evaluated
mySummary <- function (data,glm) {
  out <-
    c(Metrics::logLoss(data$actual,data$predicted),
      glm$aic)
  names(out) <- c("logLoss","AIC")
  out
}

# exceedance 100 Bq/m³ ->f(GRP)
emp.prob_GRP100 <- matrix(ncol=7,nrow=100)                   # empirical exceedance frequency
mean.prob_GRP100 <- matrix(ncol=7,nrow=100)                  # predicted probability
in.bin_GRP100 <- matrix(ncol=7,nrow=100)                     # observations in bin
in.bin.exceed_GRP100 <- matrix(ncol=7,nrow=100)
metrics_GRP100<-matrix(ncol=2,nrow=100)
#exceedance 100 Bq/m³; predictor GRP
for (i in 1:100){                                 # i=100 20*5fold cross-validation
  # fit glm model to training data
  glm <- glm(Exc_100~GRP,data=IRC[sb100[[i]],],            
             family=binomial(link="logit"))
  
  #calculate exceedance probability for test data
  predictor.data <- as.data.frame(IRC[-sb100[[i]],"GRP"])  # extract indices for test data
  colnames(predictor.data) <- "GRP"                        # name the predictor
  test.prob <- predict(glm,newdata=predictor.data,type="response")
  test.data <- as.data.frame(cbind(as.numeric(IRC[-sb100[[i]],"Exc_100"])-1,as.vector(test.prob)))
  colnames(test.data)<-c("actual","predicted")
  
  # calculate summary metrics
  metrics_GRP100[i,]<-mySummary(test.data,glm)
  
  #calculate empirical exceedance probability for all bins
  for (j in 1:7){
    indices <- which(test.data$predicted >= bin.thresh100[j] & 
                       test.data$predicted < bin.thresh100[j+1])  # find indices of test observations in bin
    in.bin_GRP100[i,j] <- length(test.data$predicted[indices])                # number of test observations in bin
    in.bin.exceed_GRP100[i,j] <- length(which(test.data$actual[indices]==1))  # number of test observations in bin exceeding 100/300 Bq/m³
    mean.prob_GRP100[i,j] <- mean(test.data$predicted[indices])                    # mean predicted probability in bin
    emp.prob_GRP100[i,j] <- in.bin.exceed_GRP100[i,j]/in.bin_GRP100[i,j]                    # mean predicted probability in bin
  }
}

# exceedance 300 Bq/m³ ->f(GRP)
emp.prob_GRP300 <- matrix(ncol=7,nrow=100)                   # empirical exceedance frequency
mean.prob_GRP300 <- matrix(ncol=7,nrow=100)                  # predicted probability
in.bin_GRP300 <- matrix(ncol=7,nrow=100)                     # observations in bin
in.bin.exceed_GRP300 <- matrix(ncol=7,nrow=100)
metrics_GRP300<-matrix(ncol=2,nrow=100)
#exceedance 100 Bq/m³; predictor GRP
for (i in 1:100){                                 # i=100 20*5fold cross-validation
  # fit glm model to training data
  glm <- glm(Exc_300~GRP,data=IRC[sb300[[i]],],            
             family=binomial(link="logit"))
  
  #calculate exceedance probability for test data
  predictor.data <- as.data.frame(IRC[-sb300[[i]],"GRP"])  # extract indices for test data
  colnames(predictor.data) <- "GRP"                        # name the predictor
  test.prob <- predict(glm,newdata=predictor.data,type="response")
  test.data <- as.data.frame(cbind(as.numeric(IRC[-sb300[[i]],"Exc_300"])-1,as.vector(test.prob)))
  colnames(test.data)<-c("actual","predicted")
  
  # calculate summary metrics
  metrics_GRP300[i,]<-mySummary(test.data,glm)
  
  #calculate empirical exceedance probability for all bins
  for (j in 1:7){
    indices <- which(test.data$predicted >= bin.thresh300[j] & 
                       test.data$predicted < bin.thresh300[j+1])  # find indices of test observations in bin
    in.bin_GRP300[i,j] <- length(test.data$predicted[indices])                # number of test observations in bin
    in.bin.exceed_GRP300[i,j] <- length(which(test.data$actual[indices]==1))  # number of test observations in bin exceeding 100/300 Bq/m³
    mean.prob_GRP300[i,j] <- mean(test.data$predicted[indices])                    # mean predicted probability in bin
    emp.prob_GRP300[i,j] <- in.bin.exceed_GRP300[i,j]/in.bin_GRP300[i,j]                    # mean predicted probability in bin
  }
}

# exceedance 100 Bq/m³ ->f(GRHI)
emp.prob_GRHI100 <- matrix(ncol=7,nrow=100)                   # empirical exceedance frequency
mean.prob_GRHI100 <- matrix(ncol=7,nrow=100)                  # predicted probability
in.bin_GRHI100 <- matrix(ncol=7,nrow=100)                     # observations in bin
in.bin.exceed_GRHI100 <- matrix(ncol=7,nrow=100)
metrics_GRHI100<-matrix(ncol=2,nrow=100)
#exceedance 100 Bq/m³; predictor GRHI
for (i in 1:100){                                 # i=100 20*5fold cross-validation
  # fit glm model to training data
  glm <- glm(Exc_100~GRHI,data=IRC[sb100[[i]],],            
             family=binomial(link="logit"))
  
  #calculate exceedance probability for test data
  predictor.data <- as.data.frame(IRC[-sb100[[i]],"GRHI"])  # extract indices for test data
  colnames(predictor.data) <- "GRHI"                        # name the predictor
  test.prob <- predict(glm,newdata=predictor.data,type="response")
  test.data <- as.data.frame(cbind(as.numeric(IRC[-sb100[[i]],"Exc_100"])-1,as.vector(test.prob)))
  colnames(test.data)<-c("actual","predicted")
  
  # calculate summary metrics
  metrics_GRHI100[i,]<-mySummary(test.data,glm)
  
  #calculate empirical exceedance probability for all bins
  for (j in 1:7){
    indices <- which(test.data$predicted >= bin.thresh100[j] & 
                       test.data$predicted < bin.thresh100[j+1])  # find indices of test observations in bin
    in.bin_GRHI100[i,j] <- length(test.data$predicted[indices])                # number of test observations in bin
    in.bin.exceed_GRHI100[i,j] <- length(which(test.data$actual[indices]==1))  # number of test observations in bin exceeding 100/100 Bq/m³
    mean.prob_GRHI100[i,j] <- mean(test.data$predicted[indices])                    # mean predicted probability in bin
    emp.prob_GRHI100[i,j] <- in.bin.exceed_GRHI100[i,j]/in.bin_GRHI100[i,j]                    # mean predicted probability in bin
  }
}
setwd("W:/UR_intern/UR2/SW1-1/Mitarbeiter/Petermann/Paper/Mapping Indoor Rn - geogenic hazard and actual risk/GitHub/Data/glm model evaluation test performance/")
# write.csv(mean.prob_GRHI100,"Prob Exc100 predicted in classes.csv")
# write.csv(emp.prob_GRHI100,"Prob Exc100 empirical in classes.csv")

# exceedance 300 Bq/m³ ->f(GRHI)
emp.prob_GRHI300 <- matrix(ncol=7,nrow=100)                   # empirical exceedance frequency
mean.prob_GRHI300 <- matrix(ncol=7,nrow=100)                  # predicted probability
in.bin_GRHI300 <- matrix(ncol=7,nrow=100)                     # observations in bin
in.bin.exceed_GRHI300 <- matrix(ncol=7,nrow=100)
metrics_GRHI300<-matrix(ncol=2,nrow=100)
#exceedance 100 Bq/m³; predictor GRHI
for (i in 1:100){                                 # i=100 20*5fold cross-validation
  # fit glm model to training data
  glm <- glm(Exc_300~GRHI,data=IRC[sb300[[i]],],            
             family=binomial(link="logit"))
  
  #calculate exceedance probability for test data
  predictor.data <- as.data.frame(IRC[-sb300[[i]],"GRHI"])  # extract indices for test data
  colnames(predictor.data) <- "GRHI"                        # name the predictor
  test.prob <- predict(glm,newdata=predictor.data,type="response")
  test.data <- as.data.frame(cbind(as.numeric(IRC[-sb300[[i]],"Exc_300"])-1,as.vector(test.prob)))
  colnames(test.data)<-c("actual","predicted")
  
  # calculate summary metrics
  metrics_GRHI300[i,]<-mySummary(test.data,glm)
  
  #calculate empirical exceedance probability for all bins
  for (j in 1:7){
    indices <- which(test.data$predicted >= bin.thresh300[j] & 
                       test.data$predicted < bin.thresh300[j+1])  # find indices of test observations in bin
    in.bin_GRHI300[i,j] <- length(test.data$predicted[indices])                # number of test observations in bin
    in.bin.exceed_GRHI300[i,j] <- length(which(test.data$actual[indices]==1))  # number of test observations in bin exceeding 100/300 Bq/m³
    mean.prob_GRHI300[i,j] <- mean(test.data$predicted[indices])                    # mean predicted probability in bin
    emp.prob_GRHI300[i,j] <- in.bin.exceed_GRHI300[i,j]/in.bin_GRHI300[i,j]                    # mean predicted probability in bin
  }
}
setwd("W:/UR_intern/UR2/SW1-1/Mitarbeiter/Petermann/Paper/Mapping Indoor Rn - geogenic hazard and actual risk/GitHub/Data/glm model evaluation test performance/")
# write.csv(mean.prob_GRHI300,"Prob Exc300 predicted in classes.csv")
# write.csv(emp.prob_GRHI300,"Prob Exc300 empirical in classes.csv")

#produce boxplots

#install.packages("ggplot2")
#install.packages("gridExtra")
library(ggplot2)
library(gridExtra)
colnames(metrics_GRP100)<-c("logLoss","AIC")
colnames(metrics_GRP300)<-c("logLoss","AIC")
colnames(metrics_GRHI100)<-c("logLoss","AIC")
colnames(metrics_GRHI300)<-c("logLoss","AIC")

logLoss100<-data.frame(logLoss=c(metrics_GRP100[,1],metrics_GRHI100[,1]),ind=c(rep("GRP",100),rep("GRHI",100)))
logLoss300<-data.frame(logLoss=c(metrics_GRP300[,1],metrics_GRHI300[,1]),ind=c(rep("GRP",100),rep("GRHI",100)))
AIC100<-data.frame(AIC=c(metrics_GRP100[,2],metrics_GRHI100[,2]),ind=c(rep("GRP",100),rep("GRHI",100)))
AIC300<-data.frame(AIC=c(metrics_GRP300[,2],metrics_GRHI300[,2]),ind=c(rep("GRP",100),rep("GRHI",100)))

logLoss100.plot <- ggplot(logLoss100, aes(x=ind, y=logLoss,fill=ind)) + 
  geom_boxplot()+
  labs(title="100 Bq/m³ exceedance",x="", y = "logLoss")+
  theme(legend.position="none")

logLoss300.plot <- ggplot(logLoss300, aes(x=ind, y=logLoss,fill=ind)) + 
  geom_boxplot()+  
  labs(title="300 Bq/m³ exceedance",x="", y = "logLoss")+
  theme(legend.position="none")

AIC100.plot <- ggplot(AIC100, aes(x=ind, y=AIC,fill=ind)) + 
  geom_boxplot()+
  labs(title="100 Bq/m³ exceedance",x="", y = "AIC")+
  theme(legend.position="none")

AIC300.plot <- ggplot(AIC300, aes(x=ind, y=AIC,fill=ind)) + 
  geom_boxplot()+
  labs(title="300 Bq/m³ exceedance",x="", y = "AIC")+
  theme(legend.position="none")

grid.arrange(logLoss100.plot,logLoss300.plot,AIC100.plot,AIC300.plot,nrow=2,ncol=2)



