#--------------------------------------------------------------------------------
#Figure S7 - Sliding window results
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
source("/Users/marcruizisagales/Documents/GitHub/Climate-baleen-plate-isotopes/analysis/4_sliding_window_nitrogen.R")
source("/Users/marcruizisagales/Documents/GitHub/Climate-baleen-plate-isotopes/analysis/6_sliding_window_carbon.R")

# 3. Sliding window results
set.seed(123)

# 4. Plot sliding window results

# nitrogen

a<-head(Bp_dN_clim_m_rel[[1]]$Dataset, 20) ## Lin_func Mean NAO, close window-open window (34-31) (beta coefficient= -0.3619122)
a <- a %>% filter(deltaAICc < (a$deltaAICc[1])+2)
nrow_a <- nrow(a)
a$Index <- c(rep("NAO", nrow_a))
b<-head(Bp_dN_clim_m_rel[[2]]$Dataset, 20) ## Lin_func Mean AMO, close window-open window (52-35) (beta coefficient= -4.377916)
b <- b %>% filter(deltaAICc < (a$deltaAICc[1])+2)
nrow_b <- nrow(b)
b$Index <- c(rep("AMO", nrow_b))

data <- rbind(a,b)

levels=c("NAO","AMO")

data$Index <- factor(data$Index, levels=rev(c("NAO","AMO")))

m1 <- ggplot(data) +
  geom_vline(aes(xintercept = data$Index), color = "gray", linetype = "dotted") +
  geom_linerange(aes(x = data$Index, ymin = WindowOpen, ymax = WindowClose), color = '#74A9CF', size = 20, alpha = 0.9) +  # Afegeix les línies de separació
  geom_hline(yintercept = seq(1,48,1), lty = 1, lwd = 0.1, colour = "white") +
  geom_hline(yintercept = 0, lty = 1, lwd = 1, colour = "black") +
  coord_flip() +
  labs(x = "", y = "Months from sample", title = "") +
  theme_article(base_size = 15) + ylim(36,0)+
  scale_x_discrete(breaks=levels,labels=c("Avg NAO index","Avg AMO index")) + standard_theme

m1

# carbon

a1<-head(Bp_dC_clim_m_rel[[2]]$Dataset, 20) ## Lin_func Mean NAO, close window-open window (34-31) (beta coefficient= -0.3619122)
a1 <- a1 %>% filter(deltaAICc < (a1$deltaAICc[1])+2)
nrow_a1 <- nrow(a1)
a1$Index <- c(rep("AMO", nrow_a1))

data1 <- rbind(a1)

levels1=c("AMO")

data1$Index <- factor(data1$Index, levels=c(levels1))

m2 <- ggplot(data1) +
  geom_vline(aes(xintercept = data1$Index), color = "gray", linetype = "dotted") +
  geom_linerange(aes(x = data1$Index, ymin = WindowOpen, ymax = WindowClose), color = '#FC8D59', size = 20, alpha =0.9) +  # Afegeix les línies de separació
  geom_hline(yintercept = seq(1,48,1), lty = 1, lwd = 0.1, colour = "white") +
  geom_hline(yintercept = 0, lty = 1, lwd = 1, colour = "black") +
  coord_flip() +
  labs(x = "", y = "Months from sample", title = "") +
  theme_article(base_size = 15) + ylim(36,0)+
  scale_x_discrete(breaks=levels1,labels=c("Avg AMO index")) + standard_theme

m1+m2

# Toghether 

m3 <- ggplot(data) +
  geom_vline(aes(xintercept = Index), color = "gray", linetype = "dotted") +
  geom_linerange(aes(x = Index, ymin = WindowOpen, ymax = WindowClose), color = '#74A9CF', size = 20, alpha = 0.9) +  # Afegeix les línies de separació
  geom_hline(yintercept = seq(1,48,1), lty = 1, lwd = 0.1, colour = "white") +
  geom_hline(yintercept = 0, lty = 1, lwd = 1, colour = "black") +
  geom_linerange(data=data1, aes(x = Index, ymin = WindowOpen, ymax = WindowClose), color = '#FC8D59', size = 20, alpha =0.9) +  # Afegeix les línies de separació
  coord_flip() +
  labs(x = "", y = "Time to baleen segment (months)", title = "") +
  theme_article(base_size = 15) + ylim(36,0)+
  scale_x_discrete(breaks=levels,labels=c("Mean NAO index","Mean AMO index")) + standard_theme
m3

ggsave("/Users/marcruizisagales/Documents/GitHub/Climate-baleen-plate-isotopes/png/Figure_S7_Sliding_window_results.png", last_plot(), dpi=300,  width = 20, height = 15, units = "cm")
