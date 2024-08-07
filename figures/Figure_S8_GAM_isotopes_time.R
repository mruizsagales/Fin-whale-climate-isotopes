#--------------------------------------------------------------------------------
# Figure S8 - Generalized Additive Models of stable isotope values trough time
#--------------------------------------------------------------------------------

# 1. Load libraries
library(tidyverse)
library(lubridate)
library(ggrepel)
library(readxl)
library(ggplot2)
library(egg)
library(zoo)
library(TSA)
library(RColorBrewer)
library(patchwork)
library(plyr)
library(dplyr)
library(lme4)
library(readr)
library(mgcv)
library(tidymv)

# 2. Import stable isotope data
df <- read_excel("~/Desktop/Doctorat/Analisis_isotops_barbes/Projecte_barbes_clima/Environmental_vars/All_merged_time_calibrated_2013_to_2022_Suess_cor.xlsx") # Import Suess-corrected dataset of stable isotope data along the baleen plate of fin whales 

# 3. Make and plot gamm models

# nitrogen
df$Whale <- as.factor(df$Whale)
df$Year_from_sample_date <- as.numeric(df$Year_from_sample_date)
d13c.cor_Year_Whale_gam <- mgcv::gam(d13c.cor ~  s(month(df$year_rev)) + s(year(df$year_rev), k=5) + s(Whale, bs = "re"), data= df, select=TRUE, method = 'REML')
summary(d13c.cor_Year_Whale_gam)

# plot nitrogen
par(mar=c(5,5,5,5))
plot.gam(d13c.cor_Year_Whale_gam)


dev.new()
png(filename="/Users/marcruizisagales/Documents/GitHub/Climate-baleen-plate-isotopes/png/Figure_S8_GAM_d13c.cor_time.png", width = 1000, height = 480)
mgcv::plot.gam(d13c.cor_Year_Whale_gam, xlab="Estimated date (years)", ylab=" s(years, 4.77)", main= expression(paste(delta^13, "C corr.")),pages = 1)
dev.off()

# carbon
df$Whale <- as.factor(df$Whale)
df$Year_from_sample_date <- as.numeric(df$Year_from_sample_date)
dN_Year_Whale_gam <- mgcv::gam(dN ~  s(month(df$year_rev)) + s(year(df$year_rev),k=5) + s(Whale, bs = "re"), data= df, select=TRUE, method = 'REML')
summary(dN_Year_Whale_gam)

# plot carbon
par(mar=c(5,5,5,5))
plot.gam(dN_Year_Whale_gam)


dev.new()
png(filename="/Users/marcruizisagales/Documents/GitHub/Climate-baleen-plate-isotopes/png/Figure_S8_GAM_dN_time.png", width = 1000, height = 480)
mgcv::plot.gam(dN_Year_Whale_gam, xlab="Estimated date (years)", ylab=" s(years, 4.09)", main= expression(paste(delta^15, "N")), pages=1)
dev.off()