#--------------------------------------------------------------------------------
#Figure S5 - Nitrogen and carbon stable isotope oscillations for each individual
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

standard_theme <- theme(
  #axis.title = element_text(size = 20),   # Axis titles
  #axis.text = element_text(size = 18),    # Axis text (tick labels)
  #plot.title = element_text(size = 16),   # Plot title
  #legend.title = element_text(size = 14), # Legend title
  legend.text = element_text(size = 12),  # Legend text
  #strip.text = element_text(size = 15),    # Facet labels, 
  legend.key.size = unit(0.5, 'cm'), 
  legend.key.height = unit(0.5, 'cm'), 
  legend.key.width = unit(0.5, 'cm'),
  aspect.ratio = 2/4,
  legend.title = element_blank(), 
  legend.position = c(0.95,0.1)
)

# 2. Import data
fw.data_d13cor <- read_excel("/Users/marcruizisagales/Documents/GitHub/Climate-baleen-plate-isotopes/data/All_merged_time_calibrated_2013_to_2022_Suess_cor.xlsx")

# 3. Plot N and C oscillations per whale

legend <- c("dC" = "#FC8D59", "dN" = "#74A9CF")

point_size=1.5
isop <- ggplot(fw.data_d13cor, aes(x = Cm)) + 
  theme_article(base_size = 20) + 
  xlab("Distance from gingiva (Cm)") + theme_article() + 
  scale_x_continuous(breaks=c(-10, 0, 10,20,30,40,50))

isop <- isop + geom_line(aes(y = dC), col= "#FC8D59") + 
  geom_point(aes(y = dC, fill= "dC"), size = point_size, shape = 21, col="white") +
  ylab(expression(paste(delta^{13}, "C (\u2030)"))) + facet_wrap(vars(Whale), ncol=8)

# add the d15N data but need to scale it to d13C
a <- mean(na.omit(fw.data_d13cor$dN))
a
b <- mean(na.omit(fw.data_d13cor$d13c.cor))
b

v_shift <- a - b

#v-shift (ind F13065) = v_shift

isop <- isop + geom_line(aes(y = dN - v_shift), col = "#74A9CF") + 
  geom_point(aes(y = (dN - v_shift), fill = "dN"), size = point_size, 
             shape = 21, col="white")+ facet_wrap(vars(Whale), ncol = 6) + 
  scale_y_continuous(sec.axis = sec_axis(~.+v_shift, name = expression(paste(delta^{15}, "N (\u2030)")))) 

its_per_whale <- isop + theme_article(base_size = 15) + scale_fill_manual(values = legend, labels = c(expression(paste(delta^{13}, "C (\u2030)")), expression(paste(delta^{15}, "N (\u2030)")))) + standard_theme
print(its_per_whale)

ggsave("/Users/marcruizisagales/Documents/GitHub/Climate-baleen-plate-isotopes/png/Figure_S5_Isotopic_oscillations_per_whale.png", last_plot(), dpi=300,  width = 30, height = 15, units = "cm")

