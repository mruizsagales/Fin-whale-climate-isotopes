##########################################################################################
### (2) Model fitting and selection R code for:

### "Timing and intensity of weather events shape nestling development strategies 
###  in three alpine breeding songbirds"

### DR de Zwaan, A Drake, JL Greenwood, and K Martin (2020)
### Frontiers in Ecology and Evolution	

### doi: 10.3389/fevo.2020.570034

### Run with R version 3.6.3
########################################################################################## 

### Set working directory

#setwd("C:/PATH")


### Load required packages

library(plyr)
library(dplyr)
library(mgcViz)
library(mgcv)
library(lme4)
library(lmerTest)
library(olsrr)
library(piecewiseSEM)
library(ggplot2)
library(tidymv)
library(doBy)
library(broom)


### Read in HOLA, DEJU, and SAVS data produced in 01_de Zwaan et al_2020_Sliding window R code

Bp_data <- read.csv("/Users/marcruizisagales/Desktop/Doctorat/Analisis_isotops_barbes/Projecte_barbes_clima/Environmental_vars/Bp_isotopic_weather_summarized_final_carbon_13_nov_mean_seasonality.csv")


###########################################################################################
### (1) Bp
###########################################################################################


######################
### Prepare data
######################

### Scale weather variables

head(Bp_data)

# Calculate the mean and standard deviation of the response variable
Bp_data_MEAN <- mean(Bp_data$d13c.cor, na.rm = TRUE)
Bp_data_SD <- sd(Bp_data$d13c.cor, na.rm = TRUE)

# Standardize the response variable
Bp_data$d13c.cor <- (Bp_data$d13c.cor - Bp_data_MEAN) / Bp_data_SD


Bp_sc <- scale(Bp_data[,c(97:111)], center=TRUE, scale=TRUE) #comprovar que s'agafen les variables

### Recombine with non-standardized variables

Bp_sc <- as.data.frame(Bp_sc)

Bp_sc$dN <- Bp_data$dN
Bp_sc$dC <- Bp_data$dC
Bp_sc$d13c.cor <- Bp_data$d13c.cor

Bp_sc$year_f <- Bp_data$year_f
Bp_sc$whale_id <- Bp_data$whale_id
Bp_sc$Species <- Bp_data$Species
Bp_sc$julian_day <- Bp_data$julian_day
Bp_sc$cohort <- Bp_data$cohort
Bp_sc$Data_capt <- Bp_data$Data_capt
Bp_sc$Length <- Bp_data$Length
Bp_sc$Status <- Bp_data$Status
Bp_sc$Sex <- Bp_data$Sex


### Convert whale id into factor

Bp_sc$whale_id <- factor(Bp_sc$whale_id)


### Convert year to a factor

Bp_sc$year_f <- factor(Bp_sc$year_f)


### Test correlation between avg temp and sub 10 hrs

#cor.test(Bp_sc$min_NAO_dN, Bp_sc$var_NAO_dN) # -0.93

### Note: correlation so high that we will continue with only average temperature, not sub 10 hrs.



##############################
### (A) Bp: dN
##############################

### Note:
# See model selection methods for details. Terms penalized out of the model are removed
# by commenting out. 
# The number behind certain variables indicates the model it was included in.

# If all terms are considered linear, then a linear model is fit below.

##############
### HOLA wing GAM
##############

Bp_dN <- mgcv::gam(d13c.cor ~  poly(julian_day,2) + 
                     s(mean_NAO_dC,k=3) +
                     s(mean_AMO_dC,k=3) +
                     s(mean_AMOC_dC_long,k=3) +
                     s(mean_AMOC_dC_short,k=3) +
                     #s(mean_SST_ice_dC,k=3) +
                     #s(mean_SST_atl_dC,k=3) +
                     #s(mean_CHL_ice_dC,k=3) +
                     #s(mean_CHL_atl_dC,k=3) +
                     s(whale_id, bs = "re"),
                   data= Bp_sc, select=TRUE, method = 'REML')

#table(Bp_sc$max_THETAO_dN_long)
# k per defecte és deu, mirar si length(table<10) i baixar la k en cas afirmatiu
#eliminar variables de edf de la mes baixa a la més alta fins a 0.7 (provar tmb amb p-valor)
summary(Bp_dN) ## 54.8% deviance explained
plot.gam(Bp_dN) #treure observació del max_AO inferior a -2
#mirar parametre layout i AiC



### Linear wing model
install.packages("lmerTest")
library(lmerTest)


Bp_dC_lm <- lmer(d13c.cor ~ poly(julian_day,2) + 
                   mean_NAO_dC +
                   mean_AMO_dC +
                   mean_AMOC_dC_long +
                   mean_AMOC_dC_short +
                   #mean_SST_ice_dC +
                   #mean_SST_atl_dC +
                   #mean_CHL_ice_dC +
                   #mean_CHL_atl_dC +
                   (1|whale_id), data=Bp_sc)

step(Bp_dC_lm)
m2 <- get_model(step(Bp_dC_lm))
summary(m2)
vif(m2)
plot(m2)
qqnorm(resid(m2)) +
  qqline(resid(m2),probs = c(0.1, 0.9))

library(nortest)
nortest::ad.test(resid(m1))
anders(resid(m1))



summary(Bp_dN_lm)

# Results
confint(Bp_dN_lm)
rsquared(Bp_dN_lm) # 56.4% R2
rsquared(m2) # 56.4% R2

install.packages("effects")
library(effects)

plot(allEffects(m2))
library(car)
vif(m2)


RColorBrewer::brewer.pal(5, "RdYlBu")

effect_plot <- effect("mean_AMO_dC", m2)
ggplot(as.data.frame(effect_plot), aes(x = mean_AMO_dC, y = fit)) +
  geom_line(color = "#FC8D59", linewidth=2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#FC8D59", alpha = 0.5,color = "black", linetype = "dotted") +
  geom_line(aes(ymin = lower, ymax = upper)) +
  ylab(expression(paste("Relative ", delta^{15}, "N (\u2030) values (SD)"))) + xlab("Avg AMO index (10-9 months)") +
  theme_classic() #+geom_point(data=Bp_sc, aes(x = Bp_sc$max_NAO_dN, y = dN))

a<-ggpredict(m2,interval = "confidence") #Confidence instead of prediction intervals can be calculated by explicitly setting interval = "confidence", in which case the random effects variances are ignored.
#a<-ggpredict(Bp_dN_lm, type = "random")
plot1 <- plot(a$mean_AMO_dC)+
  geom_line(color = "#FC8D59", linewidth=2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "#FC8D59", alpha = 0.5,color = "black", linetype = "dotted") +
  geom_line(aes(ymin = conf.low, ymax = conf.high)) +
  ylab(expression(paste("Relative ", delta^{13}, "C values (SD)"))) + xlab(expression(atop("Avg AMO index", atop(italic("(10-9 months)"), "")))) +
  theme(aspect.ratio=2/4) + ggtitle("") #+geom_point(data=Bp_sc, aes(x = Bp_sc$max_AO_dN, y = dN))

plot1 <- plot1 + theme_article(base_size = 20)

library(patchwork)
plot1 + 
  plot_layout(ncol = 3)

