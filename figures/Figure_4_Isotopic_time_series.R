#--------------------------------------------------------------------------------
#Figure 4 - Time series of nitrogen and carbon baleen stable isotopes in fin whales
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

# 2. Import data
df <- read_excel("/Users/marcruizisagales/Documents/GitHub/Climate-baleen-plate-isotopes/data/All_merged_time_calibrated_2013_to_2022_Suess_cor.xlsx") # Import Suess-corrected dataset of stable isotope data along the baleen plate of fin whales 

data <- df %>%
  pivot_longer(cols = c(dN, d13c.cor), names_to = "Class", values_to = "Isotope")

data$Class = factor(data$Class, levels=c("dN","d13c.cor")) #labels = c(expression(paste(delta^{15}, "N (\u2030)")), expression(paste(delta^{13}, "C (\u2030)"))

#dataframe transformation of the dataframe

# 3. Plot time series

nb.cols <- 29
mycolors <- colorRampPalette(brewer.pal(8, "RdYlBu"))(nb.cols) #"YlGnBu"

cyl_names <- c(
  "dN" = "delta^{15}*N",
  "d13c.cor" = "delta^{13}*C"
)


df$year_rev <- as.Date(df$year_rev)
rect <- data.frame(xmin=as.Date("2018-07-01"), xmax=as.Date("2019-06-01"), ymin=-Inf, ymax=Inf)
plot_ts <- ggplot(data=data, aes(x = as.Date(year_rev), y= Isotope, color = Whale)) +
  geom_smooth(method = "loess", se = FALSE, span = 0.3, aes(fill = Whale))  + facet_wrap(vars(Class), ncol=1, scales = "free", labeller= labeller(Class = as_labeller(cyl_names, label_parsed))) +
  xlab("Year") +
  ylab(expression(paste("Isotope ratio (\u2030)"))) +
  scale_x_date(limits = c(min(df$year_rev), max(df$year_rev)+365), date_breaks = "2 year", date_labels = "%Y") +
  geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),fill="white",alpha=1,inherit.aes = FALSE) +
  geom_vline(xintercept = rect$xmin, linetype="dashed", 
             color = "darkgrey") +
  geom_vline(xintercept = rect$xmax, linetype="dashed", 
             color = "darkgrey") + labs(fill= "Whale ID", color= "Whale ID") +
  scale_fill_manual(values = mycolors) +
  scale_color_manual(values = mycolors)+ theme_article(base_size=15) + theme(aspect.ratio = 1.5/4, legend.key.size = unit(0.5, 'cm'), 
                                                legend.key.height = unit(0.5, 'cm'), 
                                                legend.key.width = unit(0.5, 'cm'),legend.title = element_text(size=10), 
                                                legend.text = element_text(size=12))

plot_ts

ggsave("/Users/marcruizisagales/Documents/GitHub/Climate-baleen-plate-isotopes/png/Figure_4_Isotopic_time_series.png", plot_ts, 
       device = png(width = 1200, height = 450))




# 4. Plot boxplots

ggplot(data, aes(as.factor(Year_from_sample_date),Isotope)) + 
  geom_boxplot(fill="#74A9CF", alpha=0.9, width=0.2) +
  geom_jitter(color="black", size=0.1, alpha=0.1, shape=1, width = 0.1) +
  facet_wrap(vars(Class), ncol=1, scales = "free", labeller= labeller(Class = as_labeller(cyl_names, label_parsed))) +  
  #stat_summary(fun.y=mean, geom="point", shape=20, size=2, color="red", fill="red") +
  xlab("Year") + 
  ylab(expression(paste("Isotope ratio (\u2030)"))) +
  theme_article(base_size = 15) + theme(aspect.ratio=2/4) 
