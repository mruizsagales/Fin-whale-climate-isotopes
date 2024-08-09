
#--------------------------------------------------------------------------------
#Figure S6 - Beta coefficients
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

standard_theme <- theme(
  #axis.title = element_text(size = 20),   # Axis titles
  #axis.text = element_text(size = 18),    # Axis text (tick labels)
  #plot.title = element_text(size = 16),   # Plot title
  #legend.title = element_text(size = 14), # Legend title
  legend.text = element_text(size = 12),  # Legend text
  #strip.text = element_text(size = 20),    # Facet labels, 
  legend.key.size = unit(0.5, 'cm'), 
  legend.key.height = unit(0.5, 'cm'), 
  legend.key.width = unit(0.5, 'cm'),
  legend.title = element_blank(),
  legend.position = c(.9, .9),
  aspect.ratio=3/4
)

# 2. Load dependencies
source("/Users/marcruizisagales/Documents/GitHub/Climate-baleen-plate-isotopes/analysis/5_statistical_analysis_nitrogen.R")
source("/Users/marcruizisagales/Documents/GitHub/Climate-baleen-plate-isotopes/analysis/7_statistical_analysis_carbon.R")

# 3. Plot Beta coefficients
set.seed(123)

dat <- data.frame(x = rnorm(100), z = rnorm(100), y1 = rnorm(100), y2 = rnorm(100), y3 = rnorm(100))

mod1 <- m1 # nitrogen model
mod2 <- m2 # carbon model


## Create data frame of model coefficients and standard errors
ce = function(model.obj) {
  extract = summary(get(model.obj))$coefficients[ ,1:2]
  return(data.frame(extract, vars=row.names(extract), model=model.obj))
}

# Run function on the three models and bind into single data frame
coefs = do.call(rbind, sapply(paste0("mod",1:2), ce, simplify=FALSE))

names(coefs)[2] = "se" 

# Using $ notation
coefs$vars

Indexs = c( "(Intercept)","poly(julian_day, 2)1", "poly(julian_day, 2)2", "NAO index","AMO index" ,"(Intercept)", "poly(julian_day, 2)1", "poly(julian_day, 2)2", "AMO index")
coefs$Indexs = Indexs

coefs<-coefs[-1,]
coefs<-coefs[-1,]
coefs<-coefs[-1,]
coefs<-coefs[-3,]
coefs<-coefs[-3,]
coefs<-coefs[-3,]
# Faceted coefficient plot
coefs_mod1 <- filter(coefs, model== "mod1")
coefs_mod2 <- filter(coefs, model== "mod2")
pd <- position_dodge(width = 0.9)
mypalette<-brewer.pal(9,"Pastel1")

coefs$Indexs <- factor(coefs$Indexs, levels = unique(coefs$Indexs))

## Reverse the order of Subclass_Name levels
coefs$Indexs <- factor(coefs$Indexs,
                                 levels=c("NAO index","AMO index"))

levels=c("NAO index","AMO index")

coefs$Indexs <- factor(coefs$Indexs, levels=rev(c("NAO index","AMO index")))


beta_coef<- ggplot2::ggplot(coefs, aes(Indexs, Estimate, group=model)) + 
  geom_hline(yintercept=0, lty=2, lwd=1, colour="lightgrey") +
  geom_errorbar(aes(ymin=Estimate - se, ymax=Estimate + se), width=0.2, position = position_dodge(width=0.5)) +
  geom_point(aes(shape=model, fill=model), position = position_dodge(0.5),size=3, color='black') +
  coord_flip() +
  scale_x_discrete(breaks=levels,labels=c("Mean NAO index","Mean AMO index")) +
  scale_shape_manual(values = c(21:22),labels = c(expression(paste(delta^15, "N")), expression(paste(delta^13, "C")))) +
  scale_fill_manual(values=c('#FC8D59','#74A9CF'),labels = c(expression(paste(delta^15, "N")), expression(paste(delta^13, "C")))) +
  guides(colour=FALSE) +
  labs(x="", y="Standardized β coefficients") +
  theme_article(base_size = 15) + 
  standard_theme

print(beta_coef)

ggsave("/Users/marcruizisagales/Documents/GitHub/Climate-baleen-plate-isotopes/png/Figure_S6_beta_coefficients.png", last_plot(), dpi=300,  width = 20, height = 15, units = "cm")
