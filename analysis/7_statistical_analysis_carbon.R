#--------------------------------------------------------------------------------
#Statistical analysis (carbon)
#--------------------------------------------------------------------------------

# adapted from DR de Zwaan, A Drake, JL Greenwood, and K Martin (2020).
# readapted by Marc Ruiz-Sagalés May 2022.

#--------------------------------------------------------------------------------

# 1. Load libraries
library(lubridate)
library(plyr)
library(dplyr)
library(lme4)
library(climwin)
library(readxl)
library(readr)
library(ggplot2)
library(tidyverse)
library(plyr)
library(mgcv)
library(lmerTest)
library(car)
library(nortest)
library(effects)
library(nlme)
library(lvmisc)
library(egg)
library(RColorBrewer)
library(ggtext)
library(ggeffects)
library(patchwork)



# 2. Import sliding window data
Bp_data <- read.csv("/Users/marcruizisagales/Documents/GitHub/Climate-baleen-plate-isotopes/data/Bp_isotopic_climate_summarised_final_carbon_mean_seasonality_for_Github.csv")

# 3.  Prepare data

head(Bp_data)

Bp_data_MEAN <- mean(Bp_data$d13c.cor, na.rm = TRUE) # Calculate the mean of the response variable
Bp_data_SD <- sd(Bp_data$d13c.cor, na.rm = TRUE) # Calculate the SD of the response variable
Bp_data$d13c.cor <- (Bp_data$d13c.cor - Bp_data_MEAN) / Bp_data_SD # Standardize the response variable


Bp_sc <- scale(Bp_data[,c(66:71)], center=TRUE, scale=TRUE) # Scale climate variables

# Recombine with non-standardized variables
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
Bp_sc$whale_id <- factor(Bp_sc$whale_id) # Convert whale id into factor
Bp_sc$year_f <- factor(Bp_sc$year_f) # Convert year to a factor


### Test correlation between climatic variables

cor.test(Bp_sc$mean_NAO_dC, Bp_sc$mean_AMO_dC) #-0.2722169 
cor.test(Bp_sc$mean_NAO_dC, Bp_sc$mean_AMOC_dC) #0.09885563
cor.test(Bp_sc$mean_AMO_dC, Bp_sc$mean_AMOC_dC) #0.1940865


# gamm model
Bp_dC_gam <- mgcv::gam(d13c.cor ~  poly(Bp_sc$julian_day,2) + 
                         #s(mean_NAO_dC,k=3) +
                         s(mean_AMO_dC,k=3) +
                         #s(mean_AMOC_dC,k=3) + 
                         s(whale_id, bs = "re"),
                       data= Bp_sc, select=TRUE, method = 'REML')

summary(Bp_dC_gam) ## 47.8% deviance explained
plot.gam(Bp_dC_gam) # plot gamm


### lmm model

Bp_dC_lm <- lmer(d13c.cor ~ 1 + poly(julian_day,2) + 
                   mean_NAO_dC +
                   mean_AMO_dC +
                   mean_AMOC_dC +
                   (1|whale_id), data=Bp_sc)

summary(Bp_dC_lm) # NAO and AMO indices are significative

# Model selection
step(Bp_dC_lm)
m2 <- get_model(step(Bp_dC_lm))
summary(m2) # AMOC index out

# Model validation
plot_model(m2) # diagnostics
cv_m2 <- loo_cv(m2, Bp_sc, "whale_id", keep = "used") # cross-validation
accuracy(m2) # accuracy
compare_accuracy(m2, cv_m2) # compare accuracy
confint(m2) # confint

# Model validation with boostraping
qqPlot(resid(m2)) # normality of the residuals
qqPlot(c(ranef(m2)$whale_id$`(Intercept)`)) # normality of the random effects
plot(m2) # variance homogeneity

f = function(m) {res= fixef(m)[4]; names(res) = names(fixef(m)[4]); res}
boost <- lme4::bootMer(m2, f, nsim=200, type="semiparametric", use.u=TRUE)
plot(boost, index= 1)

# plots

effect_plot <- effect("mean_AMO_dC", m2)
ggplot(as.data.frame(effect_plot), aes(x = mean_AMO_dC, y = fit)) +
  geom_line(color = "#FC8D59", linewidth=2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#FC8D59", alpha = 0.5,color = "black", linetype = "dotted") +
  geom_line(aes(ymin = lower, ymax = upper)) +
  ylab(expression(paste("Relative ", delta^{15}, "N (\u2030) values (SD)"))) + xlab("Avg AMO index (10-9 months)") +
  theme_classic() #+geom_point(data=Bp_sc, aes(x = Bp_sc$max_NAO_dN, y = dN))

a<-ggpredict(m2,interval = "confidence") #Confidence instead of prediction intervals can be calculated by explicitly setting interval = "confidence", in which case the random effects variances are ignored.
#a<-ggpredict(Bp_dN_lm, type = "random")
plot2 <- plot(a$mean_AMO_dC)+
  geom_line(color = "#FC8D59", linewidth=2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "#FC8D59", alpha = 0.5,color = "black", linetype = "dotted") +
  geom_line(aes(ymin = conf.low, ymax = conf.high)) +
  ylab(expression(paste("Relative ", delta^{13}, "C values (SD)"))) + xlab(expression(atop("Avg AMO index", atop(italic("(10-9 months)"), "")))) +
  theme(aspect.ratio=2/4) + ggtitle("") #+geom_point(data=Bp_sc, aes(x = Bp_sc$max_AO_dN, y = dN))

plot2 <- plot2 + theme_article(base_size = 20) + theme(aspect.ratio = 1) +
  plot_layout(ncol = 3)
plot2

ggsave("/Users/marcruizisagales/Documents/GitHub/Climate-baleen-plate-isotopes/png/Figure_7_Predicted_trends_from_lmm_carbon.png", plot2, 
       device = png(width = 1200, height = 450))
