# install.packages("raster")
# install.packages("Smisc")
# install.packages("tibble")
# install.packages("ggplot2")
# install.packages("gridExtra")
# install.packages("dplyr")
library(raster)
library(Smisc)
library(tibble)
library(ggplot2)
library(gridExtra)
library(dplyr)

##########      #########
###  load data  ###  
#########       #########

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


## GRHI data to predict at
data100 <- with(glm100.GRHI$model, tibble(GRHI = GRHI))
data300 <- with(glm300.GRHI$model, tibble(GRHI = GRHI))

## calculate confidence interval -> add fit and se.fit on the **link** scale
data100 <- bind_cols(data100, setNames(as_tibble(predict(glm100.GRHI, data100, se.fit = TRUE)[1:2]),
                                       c('fit_link','se_link')))
data300 <- bind_cols(data300, setNames(as_tibble(predict(glm300.GRHI, data300, se.fit = TRUE)[1:2]),
                                       c('fit_link','se_link')))
## create the interval and backtransform
ilink <- family(glm100.GRHI)$linkinv        #retrieve inverse link function
data100 <- mutate(data100,
                  fit_resp  = ilink(fit_link),
                  right_upr = ilink(fit_link + (1.96 * se_link)), #95 % Confidence interval
                  right_lwr = ilink(fit_link - (1.96 * se_link)))
data300 <- mutate(data300,
                  fit_resp  = ilink(fit_link),
                  right_upr = ilink(fit_link + (1.96 * se_link)), #95 % Confidence interval
                  right_lwr = ilink(fit_link - (1.96 * se_link)))

##########      #########
###  Plot Fig 4       ###  
#########       #########

GLM100_CI <- ggplot(data100, aes(x = GRHI, y = fit_resp)) +
  ggtitle("Exceedance 100 Bq/m³")+
  geom_line(lwd=2) +
  geom_rug(aes(y = glm100.GRHI$y, colour = glm100.GRHI$model$Exc_100 )) +
  scale_colour_discrete(name = 'Exceedance') +
  scale_y_log10(name="Probability",limits=c(0.001,.7),
                breaks=c(0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5))+
  theme(
    title = element_text(size = 18),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14))+
  labs(x = 'Geogenic Radon Hazard Index', y = 'Probability (IRC>100 Bq/m³)')+
  geom_ribbon(data = data100,
                aes(ymin = right_lwr, ymax = right_upr),
                alpha = 0.1)

GLM300_CI  <- ggplot(data300, aes(x = GRHI, y = fit_resp)) +
  ggtitle("Exceedance 300 Bq/m³")+
  geom_line(lwd=2) +
  geom_rug(aes(y = glm300.GRHI$y, colour = glm300.GRHI$model$Exc_300 )) +
  scale_colour_discrete(name = 'Exceedance') +
  scale_y_log10(name="Probability",limits=c(0.001,.7),
                breaks=c(0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5))+
  theme(
    title = element_text(size = 18),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14))+
  labs(x = 'Geogenic Radon Hazard Index', y = 'Probability (IRC>300 Bq/m³)')+ 
  geom_ribbon(data = data300,aes(ymin = right_lwr, ymax = right_upr),
              alpha = 0.1)

grid.arrange(GLM100_CI,GLM300_CI,ncol=2,nrow=1)
