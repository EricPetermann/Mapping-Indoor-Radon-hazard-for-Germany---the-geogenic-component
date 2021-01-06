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

library(ggplot2)
library(mgcv)

## 100 Bq/m³
#calculate probabilities of logistic regression model
fitted.response100<-predict(glm100.GRHI,type="response")
#fit ordinary GAM for smoothing
cal.mod.GRHI_100<-gam(IRC$Exc_100~s(fitted.response100,bs="ts",m=3),family="binomial") #s() defines smooths in GAM formula; bs="ts" thin plate regression splines defines smoothing basis
#apply calibration model
calibrated.GRHI_100<-predict(cal.mod.GRHI_100,type="response")
pdata.cal<-cbind.data.frame(cal=calibrated.GRHI_100,IRC)
#plot
GLM.calib100<-ggplot(pdata.cal,aes(cal,jitter(as.numeric(Exc_100)-1,0.1)))+
  scale_x_continuous("GLM: GAM corrected values",limits=c(0,1))+
  scale_y_continuous("observed occurences",limits=c(0,1))+
  geom_point()+
  geom_abline(slope=1,intercept=0)+
  geom_smooth(col="grey60",formula=y~s(x,bs="ts",m=3))


## 300 Bq/m³
fitted.response300<-predict(glm300.GRHI,type="response")
#fit ordinary GAM for smoothing
cal.mod.GRHI_300<-gam(IRC$Exc_300~s(fitted.response300,bs="ts",m=3),family="binomial") #s() defines smooths in GAM formula; bs="ts" thin plate regression splines defines smoothing basis
#apply calibration model
calibrated.GRHI_300<-predict(cal.mod.GRHI_300,type="response")
pdata.cal<-cbind.data.frame(cal=calibrated.GRHI_300,IRC)
#plot
GLM.calib300<-ggplot(pdata.cal,aes(cal,jitter(as.numeric(Exc_300)-1,0.1)))+
  scale_x_continuous("GLM: GAM corrected values",limits=c(0,1))+
  scale_y_continuous("observed occurences",limits=c(0,1))+
  geom_point()+
  geom_abline(slope=1,intercept=0)+
  geom_smooth(col="grey60",formula=y~s(x,bs="ts",m=3))


## calibration of hold-out-data (20 times repeated 5-fold spatial cross-validation)

#loading spatial block data: 20 times repeated 5-fold spatially blocked data->100 spatial blocks with 40 km block size
sb100<-loadObject("spatialBlocks IRC Exc100.R")
sb300<-loadObject("spatialBlocks IRC Exc300.R")

# defining metrics to be evaluated
mySummary <- function (data,glm) {
  out <-
    c(Metrics::logLoss(data$actual,data$predicted),
      glm$aic)
  names(out) <- c("logLoss","AIC")
  out
}

# define ranges for bins
bin.thresh300 <- c(0,0.02,0.04,0.06,0.08,0.1,0.15,1) # 300 Bq/m³ exceedance
bin.thresh100 <- c(0,0.1,0.15,0.2,0.25,0.3,0.45,1)   # 100 Bq/m³ exceedance

# exceedance 100 Bq/m³ ->f(GRHI)
emp.prob_GRHI100.cal <- matrix(ncol=7,nrow=100)                   # empirical exceedance frequency
mean.prob_GRHI100.cal <- matrix(ncol=7,nrow=100)                  # predicted probability
in.bin_GRHI100.cal <- matrix(ncol=7,nrow=100)                     # observations in bin
in.bin.exceed_GRHI100.cal <- matrix(ncol=7,nrow=100)
metrics_GRHI100.cal <- matrix(ncol=2,nrow=100)
#exceedance 100 Bq/m³; predictor GRHI
for (i in 1:100){                                 # i=100 20*5fold cross-validation
  # fit glm model to training data
  glm <- glm(Exc_100~GRHI,data=IRC[sb100[[i]],],            
             family=binomial(link="logit"))
  
  #calculate exceedance probability for test data
  predictor.data <- as.data.frame(IRC[-sb100[[i]],"GRHI"])  # extract indices for test data
  colnames(predictor.data) <- "GRHI"                        # name the predictor
  test.prob <- as.data.frame(predict(glm,newdata=predictor.data,type="response"))
  colnames(test.prob)="fitted.response100"
  test.prob.cal<-predict(cal.mod.GRHI_100,newdata=test.prob,type="response")
  test.data <- as.data.frame(cbind(as.numeric(IRC[-sb100[[i]],"Exc_100"])-1,as.vector(test.prob.cal)))
  colnames(test.data)<-c("actual","predicted")
  
  # calculate summary metrics
  metrics_GRHI100.cal[i,]<-mySummary(test.data,glm)
  
  #calculate empirical exceedance probability for all bins
  for (j in 1:7){
    indices <- which(test.data$predicted >= bin.thresh100[j] & 
                       test.data$predicted < bin.thresh100[j+1])               # find indices of test observations in bin
    in.bin_GRHI100.cal[i,j] <- length(test.data$predicted[indices])                # number of test observations in bin
    in.bin.exceed_GRHI100.cal[i,j] <- length(which(test.data$actual[indices]==1))  # number of test observations in bin exceeding 100/100 Bq/m³
    mean.prob_GRHI100.cal[i,j] <- mean(test.data$predicted[indices])               # mean predicted probability in bin
    emp.prob_GRHI100.cal[i,j] <- in.bin.exceed_GRHI100.cal[i,j]/in.bin_GRHI100.cal[i,j]    # mean predicted probability in bin
  }
}
setwd("W:/UR_intern/UR2/SW1-1/Mitarbeiter/Petermann/Paper/Mapping Indoor Rn - geogenic hazard and actual risk/GitHub/Data/glm model evaluation test performance/")
# write.csv(mean.prob_GRHI100.cal,"Prob Exc100 predicted in classes_calibrated.csv")
# write.csv(emp.prob_GRHI100.cal,"Prob Exc100 empirical in classes_calibrated.csv")

# exceedance 300 Bq/m³ ->f(GRHI)
emp.prob_GRHI300.cal <- matrix(ncol=7,nrow=100)                   # empirical exceedance frequency
mean.prob_GRHI300.cal <- matrix(ncol=7,nrow=100)                  # predicted probability
in.bin_GRHI300.cal <- matrix(ncol=7,nrow=100)                     # observations in bin
in.bin.exceed_GRHI300.cal <- matrix(ncol=7,nrow=100)
metrics_GRHI300.cal <-matrix(ncol=2,nrow=100)
#exceedance 100 Bq/m³; predictor GRHI
for (i in 1:100){                                 # i=100 20*5fold cross-validation
  # fit glm model to training data
  glm <- glm(Exc_300~GRHI,data=IRC[sb300[[i]],],            
             family=binomial(link="logit"))
  
  #calculate exceedance probability for test data
  predictor.data <- as.data.frame(IRC[-sb300[[i]],"GRHI"])  # extract indices for test data
  colnames(predictor.data) <- "GRHI"                        # name the predictor
  test.prob <- as.data.frame(predict(glm,newdata=predictor.data,type="response"))
  colnames(test.prob)="fitted.response300"
  test.prob.cal<-predict(cal.mod.GRHI_300,newdata=test.prob,type="response")
  test.data <- as.data.frame(cbind(as.numeric(IRC[-sb300[[i]],"Exc_300"])-1,as.vector(test.prob.cal)))
  colnames(test.data)<-c("actual","predicted")
  
  # calculate summary metrics
  metrics_GRHI300.cal[i,]<-mySummary(test.data,glm)
  
  #calculate empirical exceedance probability for all bins
  for (j in 1:7){
    indices <- which(test.data$predicted >= bin.thresh300[j] & 
                       test.data$predicted < bin.thresh300[j+1])  # find indices of test observations in bin
    in.bin_GRHI300.cal[i,j] <- length(test.data$predicted[indices])                # number of test observations in bin
    in.bin.exceed_GRHI300.cal[i,j] <- length(which(test.data$actual[indices]==1))  # number of test observations in bin exceeding 100/300 Bq/m³
    mean.prob_GRHI300.cal[i,j] <- mean(test.data$predicted[indices])                    # mean predicted probability in bin
    emp.prob_GRHI300.cal[i,j] <- in.bin.exceed_GRHI300.cal[i,j]/in.bin_GRHI300.cal[i,j]                    # mean predicted probability in bin
  }
}
setwd("W:/UR_intern/UR2/SW1-1/Mitarbeiter/Petermann/Paper/Mapping Indoor Rn - geogenic hazard and actual risk/GitHub/Data/glm model evaluation test performance/")
# write.csv(mean.prob_GRHI300.cal,"Prob Exc300 predicted in classes_calibrated.csv")
# write.csv(emp.prob_GRHI300.cal,"Prob Exc300 empirical in classes_calibrated.csv")


setwd("W:/UR_intern/UR2/SW1-1/Mitarbeiter/Petermann/Paper/Mapping Indoor Rn - geogenic hazard and actual risk/GitHub/Data/glm model evaluation test performance/")
mean.prob100_predicted <- read.csv("Prob Exc100 predicted in classes.csv")[,2:8]
emp.prob100_predicted <- read.csv("Prob Exc100 empirical in classes.csv")[,2:8]
mean.prob300_predicted <- read.csv("Prob Exc300 predicted in classes.csv")[,2:8]
emp.prob300_predicted <- read.csv("Prob Exc300 empirical in classes.csv")[,2:8]
mean.prob100_calibrated <- read.csv("Prob Exc100 predicted in classes_calibrated.csv")[,2:8]
emp.prob100_calibrated <- read.csv("Prob Exc100 empirical in classes_calibrated.csv")[,2:8]
mean.prob300_calibrated <- read.csv("Prob Exc300 predicted in classes_calibrated.csv")[,2:8]
emp.prob300_calibrated <- read.csv("Prob Exc300 empirical in classes_calibrated.csv")[,2:8]


## prepare data for plots
#raw (uncalibrated data)
#install.packages ("matrixStats")
library(matrixStats)
binned_data100<-data.frame(obs=colMeans(emp.prob100_predicted,na.rm=TRUE),
                           pred= colMeans(mean.prob100_predicted,na.rm=TRUE),
                           obs90=colQuantiles(as.matrix(emp.prob100_predicted),na.rm=TRUE,prob=0.9),
                           obs10=colQuantiles(as.matrix(emp.prob100_predicted),na.rm=TRUE,prob=0.1),
                           pred90=colQuantiles(as.matrix(mean.prob100_predicted),na.rm=TRUE,prob=0.9),
                           pred10=colQuantiles(as.matrix(mean.prob100_predicted),na.rm=TRUE,prob=0.1))
binned_data300<-data.frame(obs=colMeans(emp.prob300_predicted,na.rm=TRUE),
                           pred= colMeans(mean.prob300_predicted,na.rm=TRUE),
                           obs90=colQuantiles(as.matrix(emp.prob300_predicted),na.rm=TRUE,prob=0.9),
                           obs10=colQuantiles(as.matrix(emp.prob300_predicted),na.rm=TRUE,prob=0.1),
                           pred90=colQuantiles(as.matrix(mean.prob300_predicted),na.rm=TRUE,prob=0.9),
                           pred10=colQuantiles(as.matrix(mean.prob300_predicted),na.rm=TRUE,prob=0.1))

binned_data100_cal<-data.frame(obs=colMeans(emp.prob100_calibrated,na.rm=TRUE),
                           pred= colMeans(mean.prob100_calibrated,na.rm=TRUE),
                           obs90=colQuantiles(as.matrix(emp.prob100_calibrated),na.rm=TRUE,prob=0.9),
                           obs10=colQuantiles(as.matrix(emp.prob100_calibrated),na.rm=TRUE,prob=0.1),
                           pred90=colQuantiles(as.matrix(mean.prob100_calibrated),na.rm=TRUE,prob=0.9),
                           pred10=colQuantiles(as.matrix(mean.prob100_calibrated),na.rm=TRUE,prob=0.1))
binned_data300_cal<-data.frame(obs=colMeans(emp.prob300_calibrated,na.rm=TRUE),
                           pred= colMeans(mean.prob300_calibrated,na.rm=TRUE),
                           obs90=colQuantiles(as.matrix(emp.prob300_calibrated),na.rm=TRUE,prob=0.9),
                           obs10=colQuantiles(as.matrix(emp.prob300_calibrated),na.rm=TRUE,prob=0.1),
                           pred90=colQuantiles(as.matrix(mean.prob300_calibrated),na.rm=TRUE,prob=0.9),
                           pred10=colQuantiles(as.matrix(mean.prob300_calibrated),na.rm=TRUE,prob=0.1))


#prob100
ilink <- family(cal.mod.GRHI_100)$linkinv  #get link function
data100 <- with(cal.mod.GRHI_100, tibble(fitted.response100 = fitted.response100))
## add fit and se.fit on the **link** scale
data100 <- bind_cols(data100, setNames(as_tibble(predict(cal.mod.GRHI_100, data100, se.fit = TRUE)[1:2]),
                                         c('fit_link','se_link')))
## create the interval and backtransform
data100 <- mutate(data100,
                   fit_cal  = ilink(fit_link),
                   cal_upr = ilink(fit_link + (1.282 * se_link)), #z-distribution: 80 % confidence interval (10-90 %ile)
                   cal_lwr = ilink(fit_link - (1.282 * se_link)))

#prob300
data300 <- with(cal.mod.GRHI_300, tibble(fitted.response300 = fitted.response300))
## add fit and se.fit on the **link** scale
data300 <- bind_cols(data300, setNames(as_tibble(predict(cal.mod.GRHI_300, data300, se.fit = TRUE)[1:2]),
                                       c('fit_link','se_link')))
## create the interval and backtransform
data300 <- mutate(data300,
                  fit_cal  = ilink(fit_link),
                  cal_upr = ilink(fit_link + (1.282 * se_link)), #z-Verteilungstablle: 80 % Ii (10-90)
                  cal_lwr = ilink(fit_link - (1.282 * se_link)))


pred100plot<-ggplot(data100, aes(x = fitted.response100, y = fit_cal)) +
  ggtitle("Exceedance 100 Bq/m³")+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_line(lwd=2,colour="gray")+
  geom_point(data=binned_data100,mapping=aes(x=pred,y=obs),size=3, inherit.aes = FALSE)+
  geom_errorbar(data=binned_data100,mapping=aes(x=pred,ymin=obs10, ymax=obs90), width=.05,
                position=position_dodge(.9), inherit.aes = FALSE) +
  geom_errorbarh(data=binned_data100,mapping=aes(y=obs,xmin=pred10, xmax=pred90), height=.05,
                 position=position_dodge(.9), inherit.aes = FALSE) +
  scale_y_log10(name="Observed frequency",limits=c(0.03,1),breaks=c(.01,.05,0.1,.2,.5))+
  scale_x_log10(name="Raw probability",limits=c(0.03,1),breaks=c(.01,.05,0.1,.2,.5))+
  theme_gray()+
  geom_ribbon(data = data100,
              aes(ymin = cal_lwr, ymax = cal_upr),
              alpha = 0.1)+
  geom_abline(slope=1,intercept=0)+ 
  theme(axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14))


pred300plot<-ggplot(data300, aes(x = fitted.response300, y = fit_cal)) +
  ggtitle("Exceedance 300 Bq/m³")+
  geom_line(lwd=2,colour="gray")+
  geom_point(data=binned_data300,mapping=aes(x=pred,y=obs),size=3)+
  geom_errorbar(data=binned_data300,mapping=aes(x=pred,ymin=obs10, ymax=obs90), width=.05,
                position=position_dodge(.9), inherit.aes = FALSE) +
  geom_errorbarh(data=binned_data300,mapping=aes(y=obs,xmin=pred10, xmax=pred90), height=.05,
                 position=position_dodge(.9), inherit.aes = FALSE) +
  scale_y_log10(name="Observed frequency",limits=c(0.005,.55),breaks=c(0.01,0.02,.05,0.1,.2,.5))+
  scale_x_log10(name="Raw probability",limits=c(0.005,.55),breaks=c(.01,.02,.05,0.1,.2,.5))+
  theme_gray()+
  geom_ribbon(data = data300,
              aes(ymin = cal_lwr, ymax = cal_upr),
              alpha = 0.1)+
  geom_abline(slope=1,intercept=0)+
  theme(axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14))


cal100GAMplot<-ggplot(binned_data100_cal) +
  ggtitle("Exceedance 100 Bq/m³")+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_point(data=binned_data100_cal,mapping=aes(x=pred,y=obs),size=3, inherit.aes = FALSE)+
  geom_errorbar(data=binned_data100_cal,mapping=aes(x=pred,ymin=obs10, ymax=obs90), width=.05,
                position=position_dodge(.9), inherit.aes = FALSE) +
  geom_errorbarh(data=binned_data100_cal,mapping=aes(y=obs,xmin=pred10, xmax=pred90), height=.05,
                 position=position_dodge(.9), inherit.aes = FALSE) +
  scale_y_log10(name="Observed frequency",limits=c(0.03,1),breaks=c(.01,.05,0.1,.2,.5))+
  scale_x_log10(name="Calibrated probability",limits=c(0.03,1),breaks=c(.01,.05,0.1,.2,.5))+
  theme_gray()+
  geom_abline(slope=1,intercept=0)+
  theme(axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14))


cal300GAMplot<-ggplot(binned_data300_cal) +
  ggtitle("Exceedance 300 Bq/m³")+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_point(data=binned_data300_cal,mapping=aes(x=pred,y=obs),size=3, inherit.aes = FALSE)+
  geom_errorbar(data=binned_data300_cal,mapping=aes(x=pred,ymin=obs10, ymax=obs90), width=.05,
                position=position_dodge(.9), inherit.aes = FALSE) +
  geom_errorbarh(data=binned_data300_cal,mapping=aes(y=obs,xmin=pred10, xmax=pred90), height=.05,
                 position=position_dodge(.9), inherit.aes = FALSE) +
  scale_y_log10(name="Observed frequency",limits=c(0.005,.55),breaks=c(0.01,0.02,.05,0.1,.2,.5))+
  scale_x_log10(name="Calibrated probability",limits=c(0.005,.55),breaks=c(.01,.02,.05,0.1,.2,.5))+
  theme_gray()+
  geom_abline(slope=1,intercept=0)+
  theme(axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14))


grid.arrange(pred100plot,pred300plot,cal100GAMplot,cal300GAMplot,nrow=2,ncol=2)
