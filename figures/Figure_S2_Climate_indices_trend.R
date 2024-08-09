#--------------------------------------------------------------------------------
#Plot North Atlantic climate indices
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
library(climwin)
library(readr)

# 2. Load monthly weather data
data_climate <- read_excel("/Users/marcruizisagales/Desktop/Papers/Paper_climate/Datasets_paper_clim/merged_climate_cycles_dated_ts.xlsx") #load climate data
data_climate<- as.data.frame(data_climate)
data_climate$Date <- as.Date(data_climate$date)
data_climate$date_use <- format(data_climate$Date, format= "%d/%m/%Y")

data_env <- read_excel("/Users/marcruizisagales/Desktop/Papers/Paper_climate/Datasets_paper_clim/Env_total_13_nov.xlsx") # load environmental data
data_env<- as.data.frame(data_env)
data_env$Date <- as.Date(data_env$Date)

data_climate<- merge(data_climate, data_env, by = "Date", all=TRUE)

#standardize climate indices

NAO_MEAN <- mean(data_climate$NAO_index, na.rm = TRUE) # Calculate the mean of the response variable
NAO_SD <- sd(data_climate$NAO_index, na.rm = TRUE) # Calculate the SD of the response variable
data_climate$NAO_index <- (data_climate$NAO_index - NAO_MEAN) / NAO_SD # Standardize the response variable

AMO_MEAN <- mean(data_climate$AMO_index, na.rm = TRUE) # Calculate the mean of the response variable
AMO_SD <- sd(data_climate$AMO_index, na.rm = TRUE) # Calculate the SD of the response variable
data_climate$AMO_index <- (data_climate$AMO_index - AMO_MEAN) / AMO_SD # Standardize the response variable

AMOC_MEAN <- mean(data_climate$AMOC_index, na.rm = TRUE) # Calculate the mean of the response variable
AMOC_SD <- sd(data_climate$AMOC_index, na.rm = TRUE) # Calculate the SD of the response variable
data_climate$AMOC_index <- (data_climate$AMOC_index - AMOC_MEAN) / AMOC_SD # Standardize the response variable

data_climate_filt <- data_climate[711:901,]
install.packages("egg")

RColorBrewer::brewer.pal(3, "Greys")
ggplot(data=data_climate_filt, aes(x=as.Date(data_climate_filt$date))) + 
  geom_line(aes(y=data_climate_filt$AMOC_index), color= "black", size=1)+ 
  geom_line(aes(y=data_climate_filt$AMO_index), color=  "#BDBDBD", size=1)+ 
  geom_line(aes(y=data_climate_filt$NAO_index), color= "#636363", size=1)+ 
  theme_article(base_size = 15) + 
  ylab("Climate index") + 
  xlab("Year") + 
  theme(aspect.ratio = 2/4) + 
  scale_x_date(date_breaks = "2 years", date_labels = "%Y")

ggsave("/Users/marcruizisagales/Desktop/Climateindexs.png", last_plot(), dpi = 300,  width = 20, height = 15, units = "cm")

#a <- data.frame(class= rep(c("AMO index", "NAO index","AMOC index"),4), date = c(1,2,3,4,5,6,7,8,9,10,11,12), number= rep(rnorm(1,3,2),12))
#a$class <- factor(a$class, levels= c("NAO index", "AMO index", "AMOC index"))
#col= c("#636363", "#BDBDBD", "black")
#ggplot(data=a, aes(x=date, y=number)) + geom_line(aes(color=class), size=1) + scale_color_manual(values =col) + theme_article()
#ggsave("/Users/marcruizisagales/Desktop/Climateindexs1.png", last_plot(), dpi = 300,  width = 20, height = 15, units = "cm")

library(tidyr)
data_climate_filt_bo <-data_climate_filt[,c("date", "AMO_index", "AMOC_index", "NAO_index")]
df_long <- data_climate_filt_bo %>%
  pivot_longer(cols = 2:4,
               names_to = "climate_index",
               values_to = "climate_value")

ggplot(data=df_long, aes(x=as.Date(df_long$date))) + 
  geom_line(aes(y=df_long$climate_value, color= climate_index), size=1)+ 
  
  theme_article(base_size = 15) + 
  ylab("Climate index") + 
  xlab("Year") + 
  theme(aspect.ratio = 2/4) + 
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") + facet_wrap(vars(climate_index))

# Print the restructured dataframe
print(df_long)

plot(data=data_climate_filt, aes(data_climate_filt$AMO_index,data_climate_filt$NAO_index)) + geom_point() + geom_smooth(method= "lm")
