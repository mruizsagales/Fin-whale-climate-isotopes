#--------------------------------------------------------------------------------
#Time-calibration of stable isotope oscillations in baleen plates from North Atlantic fin whales
#--------------------------------------------------------------------------------

# inspired by Clive Trueman's script.
# adapted by Andrew Jackson June 2017 & Natalie Cooper Oct 2017.
# readapted by Marc Ruiz-Sagalés May 2022.

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

df <- read_excel("/Users/marcruizisagales/Documents/GitHub/Climate-baleen-plate-isotopes/data/0_point_alignments_29_juny_total_-4.xlsx") # Import dataset of the one-centimetre-spaced stable isotope data along the baleen plate of fin whales

min_cm_minus_4 <- df %>% dplyr::group_by(Whale) %>% dplyr::summarise(min_cm = min(Cm)) # Make sure that all the stable isotope data start at the position -4 cm 

# 3. Time-calibration of the stable isotope data
#________________________________________________________________________________
# Samples from 2013 whales
#________________________________________________________________________________

df_2013 <- dplyr::filter(df, Year == 2013)
unique(df_2013$Whale) # (n=7; "F13065" "F13066" "F13068" "F13073" "F13076" "F13129" "F13083")

point_size <- 2
combined_df_2013 <- NULL;

for (i in 1:length(unique(df_2013$Whale))) { 
  df_i <- dplyr::filter(df, Whale == unique(df_2013$Whale)[i]) # We select one fin whale individual from 2013
  growth_rate <- 16  # Baleen plate growth rate (centimeters per year)
  last_sample <- unique(ymd(df_i$Data_capt)) # Date of last keratin sample (catch date)
  df_i <- 
    df_i %>%
    dplyr::mutate(days = Cm * 365.25 / growth_rate) %>% # Number of days between samples starting at Cm 0
    dplyr::mutate(days = days - days[1]) %>% # Number of days between samples starting at Cm -4
    dplyr::mutate(rev_days = days - days[length(days)]) %>% # Number of days starting at maximum Cm
    dplyr::mutate(sample_date = last_sample + rev_days) %>% # Estimated date
    dplyr::mutate(year = year(sample_date)) %>% # Estimated date year
    dplyr::mutate(year_rev = rev(sample_date)) # Estimated date retrospectively (real date)
  
  df_i$Year_from_sample_date <- format(as.Date(df_i$year_rev), "%Y") # year from the date estimated retrospectively
  df_i$Year_month <- format(as.Date(df_i$year_rev), "%Y/%m") # year and month from the date estimated retrospectively
  df_i$Month_from_sample_date <- format(as.Date(df_i$year_rev), "%m") # month from the date estimated retrospectively
  
  # plot top and bottom axis
  isop <- ggplot(df_i, aes(x = rev(Cm))) + 
    theme_classic(base_size = 16) + 
    xlab("Baleen sample (cm from youngest sample)") + 
    scale_x_continuous(breaks = df_i$Cm, 
                       labels = rev(df_i$Cm),
                       sec.axis = sec_axis(~., 
                                           breaks = rev(df_i$Cm),
                                           labels = rev(df_i$sample_date),
                                           name = "Date"))
  
  #isop <- ggplot(df_i, aes(x = rev(sample_date))) + xlab("Date") + scale_x_continuous(breaks = seq(min(df_i$sample_date), max(df_i$sample_date), by = "1 month"),labels = format(seq(min(rev(df_i$sample_date)), max(rev(df_i$sample_date)), by = "1 month"), "%b-%Y"))
  
  # plot carbon stable isotope values
  isop <- isop + geom_line(aes(y = dC), col = "black") + 
    geom_point(aes(y = dC), size = point_size, shape = 21, fill = "black", col = "black") +
    ylab(expression(paste(delta^{13}, "C (\u2030)")))+ theme_article()  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  a <- mean(na.omit(df_i$dN))# add the d15N data but need to scale it to d13C
  a
  b <- mean(na.omit(df_i$dC))# add the d15N data but need to scale it to d13C
  b
  v_shift <- a - b
  
  # plot nitrogen stable isotope values
  isop <- isop + geom_line(aes(y = dN - v_shift), col = "darkgrey") + 
    geom_point(aes(y = (dN - v_shift)), size = point_size, 
               shape = 21, fill = "grey", col = "darkgrey")
  
  # plot x and y axis
  isop <- isop + 
    scale_y_continuous(sec.axis = sec_axis(~.+v_shift, 
                                           name = expression(paste(delta^{15}, "N (\u2030)"))),limits = c(-30, -10)) + labs(title= unique(df_2013$Whale)[i])
  new_df <- df_i
  combined_df_2013 <- rbind(combined_df_2013, new_df)
  print(isop)
}

# plot time-calibrated stable isotope data for the individuals of 2013
brewer.pal(n = 7, name = "Reds") #palette for plot
unique(combined_df_2013$Whale)

plot_dN_2013 <- ggplot() + geom_line(data= combined_df_2013[combined_df_2013$Whale == "F13065",], aes(x= rev(sample_date), y= dN),col="#FEE5D9")+ 
  geom_point(data= combined_df_2013[combined_df_2013$Whale == "F13065",],aes(x= rev(sample_date), y= dN), size = point_size, shape = 21, fill = "#FEE5D9") +
  geom_line(data= combined_df_2013[combined_df_2013$Whale == "F13066",], aes(x= rev(sample_date), y= dN),col="#FCBBA1")+ 
  geom_point(data= combined_df_2013[combined_df_2013$Whale == "F13066",],aes(x= rev(sample_date), y= dN), size = point_size, shape = 21, fill = "#FCBBA1") + 
  geom_line(data= combined_df_2013[combined_df_2013$Whale == "F13068",], aes(x= rev(sample_date), y= dN),col="#FC9272") +
  geom_point(data= combined_df_2013[combined_df_2013$Whale == "F13068",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#FC9272") +
  geom_line(data= combined_df_2013[combined_df_2013$Whale == "F13073",], aes(x= rev(sample_date), y= dN),col="#FB6A4A")+ 
  geom_point(data= combined_df_2013[combined_df_2013$Whale == "F13073",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#FB6A4A") + 
  geom_line(data= combined_df_2013[combined_df_2013$Whale == "F13076",], aes(x= rev(sample_date), y= dN),col="#EF3B2C")+ 
  geom_point(data= combined_df_2013[combined_df_2013$Whale == "F13076",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#EF3B2C") +
  #geom_line(data= combined_df_2013[combined_df_2013$Whale == "F13083",], aes(x= rev(sample_date), y= dN),col="#CB181D")+ 
  #geom_point(data= combined_df_2013[combined_df_2013$Whale == "F13083",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#CB181D") +
  #geom_line(data= combined_df_2013[combined_df_2013$Whale == "F13129",], aes(x= rev(sample_date), y= dN),col="#99000D")+ 
  #geom_point(data= combined_df_2013[combined_df_2013$Whale == "F13129",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#99000D") +
  theme_article(base_size = 15) + 
  xlab("Year") + theme(aspect.ratio = 3/4) +
  ylab(expression(paste(delta^{15}, "N (\u2030)")))

plot_dC_2013 <-ggplot() + geom_line(data= combined_df_2013[combined_df_2013$Whale == "F13065",], aes(x= rev(sample_date), y= dC),col="#FEE5D9")+ 
  geom_point(data= combined_df_2013[combined_df_2013$Whale == "F13065",],aes(x= rev(sample_date), y= dC), size = point_size, shape = 21, fill = "#FEE5D9") +
  geom_line(data= combined_df_2013[combined_df_2013$Whale == "F13066",], aes(x= rev(sample_date), y= dC),col="#FCBBA1")+ 
  geom_point(data= combined_df_2013[combined_df_2013$Whale == "F13066",],aes(x= rev(sample_date), y= dC), size = point_size, shape = 21, fill = "#FCBBA1") + 
  geom_line(data= combined_df_2013[combined_df_2013$Whale == "F13068",], aes(x= rev(sample_date), y= dC),col="#FC9272") +
  geom_point(data= combined_df_2013[combined_df_2013$Whale == "F13068",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#FC9272") +
  geom_line(data= combined_df_2013[combined_df_2013$Whale == "F13073",], aes(x= rev(sample_date), y= dC),col="#FB6A4A")+ 
  geom_point(data= combined_df_2013[combined_df_2013$Whale == "F13073",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#FB6A4A") + 
  geom_line(data= combined_df_2013[combined_df_2013$Whale == "F13076",], aes(x= rev(sample_date), y= dC),col="#EF3B2C")+ 
  geom_point(data= combined_df_2013[combined_df_2013$Whale == "F13076",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#EF3B2C") +
  #geom_line(data= combined_df_2013[combined_df_2013$Whale == "F13083",], aes(x= rev(sample_date), y= dC),col="#CB181D")+ 
  #geom_point(data= combined_df_2013[combined_df_2013$Whale == "F13083",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#CB181D") +
  #geom_line(data= combined_df_2013[combined_df_2013$Whale == "F13129",], aes(x= rev(sample_date), y= dC),col="#99000D")+ 
  #geom_point(data= combined_df_2013[combined_df_2013$Whale == "F13129",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#99000D") +
  theme_article(base_size = 15) + 
  xlab("Year") + theme(aspect.ratio = 3/4) +
  ylab(expression(paste(delta^{13}, "C (\u2030)")))

plot_dN_2013/plot_dC_2013

# Average values for 2013 individuals

#Nitrogen stable isotope values

#merge_2013_filtered <- combined_df_2013 %>% dplyr::group_by(Cm) %>% dplyr::summarise(across(.cols = where(is.numeric), .fns = ~mean(dN, na.rm = TRUE))) # Average value per centimeter
#ggplot(merge_2013_filtered, aes(x=merge_2013_filtered$Cm,y=merge_2013_filtered$dN)) + geom_point() + geom_line() + theme_article() + xlab("Cm") + ylab(expression(paste(delta^{15}, "N (\u2030)")))  

#Carbon stable isotope values

#merge_2013_filtered <- combined_df_2013 %>% dplyr::group_by(Cm) %>% dplyr::summarise(across(.cols = where(is.numeric), .fns = ~mean(dC, na.rm = TRUE))) # Average value per centimeter
#ggplot(merge_2013_filtered, aes(x=merge_2013_filtered$Cm,y=merge_2013_filtered$dC)) + geom_point() + geom_line() + theme_article() + xlab("Cm") + ylab(expression(paste(delta^{13}, "C (\u2030)")))  

#________________________________________________________________________________
# Samples from 2015 whales
#________________________________________________________________________________

df_2015 <- dplyr::filter(df, Year == 2015)
unique(df_2015$Whale) # (n=9; "F15078" "F15097" "F15086" "F15083" "F15080" "F15079" "F15084" "F15087" "F15088")

point_size <- 2
combined_df_2015 <- NULL;

for (i in 1:length(unique(df_2015$Whale))) { 
  df_i <- dplyr::filter(df, Whale == unique(df_2015$Whale)[i]) # We select one fin whale individual from 2015
  growth_rate <- 16  # Baleen plate growth rate (centimeters per year)
  last_sample <- unique(ymd(df_i$Data_capt)) # Date of last keratin sample (catch date)
  df_i <- 
    df_i %>%
    dplyr::mutate(days = Cm * 365.25 / growth_rate) %>% # Number of days between samples starting at Cm 0
    dplyr::mutate(days = days - days[1]) %>% # Number of days between samples starting at Cm -4
    dplyr::mutate(rev_days = days - days[length(days)]) %>% # Number of days starting at maximum Cm
    dplyr::mutate(sample_date = last_sample + rev_days) %>% # Estimated date
    dplyr::mutate(year = year(sample_date)) %>% # Estimated date year
    dplyr::mutate(year_rev = rev(sample_date)) # Estimated date retrospectively (real date)
  
  df_i$Year_from_sample_date <- format(as.Date(df_i$year_rev), "%Y") # year from the date estimated retrospectively
  df_i$Year_month <- format(as.Date(df_i$year_rev), "%Y/%m") # year and month from the date estimated retrospectively
  df_i$Month_from_sample_date <- format(as.Date(df_i$year_rev), "%m") # month from the date estimated retrospectively
  
  # plot top and bottom axis
  isop <- ggplot(df_i, aes(x = rev(Cm))) + 
    theme_classic(base_size = 16) + 
    xlab("Baleen sample (cm from youngest sample)") + 
    scale_x_continuous(breaks = df_i$Cm, 
                       labels = rev(df_i$Cm),
                       sec.axis = sec_axis(~., 
                                           breaks = rev(df_i$Cm),
                                           labels = rev(df_i$sample_date),
                                           name = "Date"))
  
  #isop <- ggplot(df_i, aes(x = rev(sample_date))) + xlab("Date") + scale_x_continuous(breaks = seq(min(df_i$sample_date), max(df_i$sample_date), by = "1 month"),labels = format(seq(min(rev(df_i$sample_date)), max(rev(df_i$sample_date)), by = "1 month"), "%b-%Y"))
  
  # plot carbon stable isotope values
  isop <- isop + geom_line(aes(y = dC), col = "black") + 
    geom_point(aes(y = dC), size = point_size, shape = 21, fill = "black", col = "black") +
    ylab(expression(paste(delta^{13}, "C (\u2030)")))+ theme_article()  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  a <- mean(na.omit(df_i$dN))# add the d15N data but need to scale it to d13C
  a
  b <- mean(na.omit(df_i$dC))# add the d15N data but need to scale it to d13C
  b
  v_shift <- a - b
  
  # plot nitrogen stable isotope values
  isop <- isop + geom_line(aes(y = dN - v_shift), col = "darkgrey") + 
    geom_point(aes(y = (dN - v_shift)), size = point_size, 
               shape = 21, fill = "grey", col = "darkgrey")
  
  # plot x and y axis
  isop <- isop + 
    scale_y_continuous(sec.axis = sec_axis(~.+v_shift, 
                                           name = expression(paste(delta^{15}, "N (\u2030)"))),limits = c(-30, -10)) + labs(title= unique(df_2015$Whale)[i])
  new_df <- df_i
  combined_df_2015 <- rbind(combined_df_2015, new_df)
  print(isop)
}

# plot time-calibrated stable isotope data for the individuals of 2015
brewer.pal(n = 9, name = "Greens") #palette for plot
unique(combined_df_2015$Whale)

plot_dN_2015 <- ggplot() + geom_line(data= combined_df_2015[combined_df_2015$Whale == "F15078",], aes(x= rev(sample_date), y= dN),col="#F7FCF5")+ 
  geom_point(data= combined_df_2015[combined_df_2015$Whale == "F15078",],aes(x= rev(sample_date), y= dN), size = point_size, shape = 21, fill = "#F7FCF5") +
  geom_line(data= combined_df_2015[combined_df_2015$Whale == "F15097",], aes(x= rev(sample_date), y= dN),col="#E5F5E0")+ 
  geom_point(data= combined_df_2015[combined_df_2015$Whale == "F15097",],aes(x= rev(sample_date), y= dN), size = point_size, shape = 21, fill = "#E5F5E0") + 
  geom_line(data= combined_df_2015[combined_df_2015$Whale == "F15086",], aes(x= rev(sample_date), y= dN),col="#C7E9C0") +
  geom_point(data= combined_df_2015[combined_df_2015$Whale == "F15086",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#C7E9C0") +
  geom_line(data= combined_df_2015[combined_df_2015$Whale == "F15083",], aes(x= rev(sample_date), y= dN),col="#A1D99B")+ 
  geom_point(data= combined_df_2015[combined_df_2015$Whale == "F15083",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#A1D99B") + 
  geom_line(data= combined_df_2015[combined_df_2015$Whale == "F15080",], aes(x= rev(sample_date), y= dN),col="#74C476")+ 
  geom_point(data= combined_df_2015[combined_df_2015$Whale == "F15080",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#74C476") +
  geom_line(data= combined_df_2015[combined_df_2015$Whale == "F15079",], aes(x= rev(sample_date), y= dN),col="#41AB5D")+ 
  geom_point(data= combined_df_2015[combined_df_2015$Whale == "F15079",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#41AB5D") +
  geom_line(data= combined_df_2015[combined_df_2015$Whale == "F15084",], aes(x= rev(sample_date), y= dN),col="#238B45")+ 
  geom_point(data= combined_df_2015[combined_df_2015$Whale == "F15084",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#238B45") +
  geom_line(data= combined_df_2015[combined_df_2015$Whale == "F15087",], aes(x= rev(sample_date), y= dN),col="#006D2C")+ 
  geom_point(data= combined_df_2015[combined_df_2015$Whale == "F15087",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#006D2C") +
  geom_line(data= combined_df_2015[combined_df_2015$Whale == "F15088",], aes(x= rev(sample_date), y= dN),col="#00441B")+ 
  geom_point(data= combined_df_2015[combined_df_2015$Whale == "F15088",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#00441B") +
  theme_article(base_size = 15) + 
  xlab("Year") + theme(aspect.ratio = 3/4) +
  ylab(expression(paste(delta^{15}, "N (\u2030)")))

plot_dC_2015 <- ggplot() + geom_line(data= combined_df_2015[combined_df_2015$Whale == "F15078",], aes(x= rev(sample_date), y= dC),col="#F7FCF5")+ 
  geom_point(data= combined_df_2015[combined_df_2015$Whale == "F15078",],aes(x= rev(sample_date), y= dC), size = point_size, shape = 21, fill = "#F7FCF5") +
  geom_line(data= combined_df_2015[combined_df_2015$Whale == "F15097",], aes(x= rev(sample_date), y= dC),col="#E5F5E0")+ 
  geom_point(data= combined_df_2015[combined_df_2015$Whale == "F15097",],aes(x= rev(sample_date), y= dC), size = point_size, shape = 21, fill = "#E5F5E0") + 
  geom_line(data= combined_df_2015[combined_df_2015$Whale == "F15086",], aes(x= rev(sample_date), y= dC),col="#C7E9C0") +
  geom_point(data= combined_df_2015[combined_df_2015$Whale == "F15086",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#C7E9C0") +
  geom_line(data= combined_df_2015[combined_df_2015$Whale == "F15083",], aes(x= rev(sample_date), y= dC),col="#A1D99B")+ 
  geom_point(data= combined_df_2015[combined_df_2015$Whale == "F15083",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#A1D99B") + 
  geom_line(data= combined_df_2015[combined_df_2015$Whale == "F15080",], aes(x= rev(sample_date), y= dC),col="#74C476")+ 
  geom_point(data= combined_df_2015[combined_df_2015$Whale == "F15080",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#74C476") +
  geom_line(data= combined_df_2015[combined_df_2015$Whale == "F15079",], aes(x= rev(sample_date), y= dC),col="#41AB5D")+ 
  geom_point(data= combined_df_2015[combined_df_2015$Whale == "F15079",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#41AB5D") +
  geom_line(data= combined_df_2015[combined_df_2015$Whale == "F15084",], aes(x= rev(sample_date), y= dC),col="#238B45")+ 
  geom_point(data= combined_df_2015[combined_df_2015$Whale == "F15084",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#238B45") +
  geom_line(data= combined_df_2015[combined_df_2015$Whale == "F15087",], aes(x= rev(sample_date), y= dC),col="#006D2C")+ 
  geom_point(data= combined_df_2015[combined_df_2015$Whale == "F15087",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#006D2C") +
  geom_line(data= combined_df_2015[combined_df_2015$Whale == "F15088",], aes(x= rev(sample_date), y= dC),col="#00441B")+ 
  geom_point(data= combined_df_2015[combined_df_2015$Whale == "F15088",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#00441B") +
  theme_article(base_size = 15) + 
  xlab("Year") + theme(aspect.ratio = 3/4) +
  ylab(expression(paste(delta^{13}, "C (\u2030)")))

plot_dN_2015/plot_dC_2015

# Average values for 2015 individuals

#Nitrogen stable isotope values

#merge_2015_filtered <- combined_df_2015 %>% dplyr::group_by(Cm) %>% dplyr::summarise(across(.cols = where(is.numeric), .fns = ~mean(dN, na.rm = TRUE))) # Average value per centimeter
#ggplot(merge_2015_filtered, aes(x=merge_2015_filtered$Cm,y=merge_2015_filtered$dN)) + geom_point() + geom_line() + theme_article() + xlab("Cm") + ylab(expression(paste(delta^{15}, "N (\u2030)")))  

#Carbon stable isotope values

#merge_2015_filtered <- combined_df_2015 %>% dplyr::group_by(Cm) %>% dplyr::summarise(across(.cols = where(is.numeric), .fns = ~mean(dC, na.rm = TRUE))) # Average value per centimeter
#ggplot(merge_2015_filtered, aes(x=merge_2015_filtered$Cm,y=merge_2015_filtered$dC)) + geom_point() + geom_line() + theme_article() + xlab("Cm") + ylab(expression(paste(delta^{13}, "C (\u2030)")))  

#________________________________________________________________________________
# Samples from 2018 whales
#________________________________________________________________________________

df_2018 <- dplyr::filter(df, Year == 2018)
unique(df_2018$Whale) # (n=13; "F18044" "F18043" "F18038" "F18050" "F18051" "F18030" "F18054" "F18007" "F18003" "F18009" "F18025" "F18098" "F18022")

point_size <- 2
combined_df_2018 <- NULL;

for (i in 1:length(unique(df_2018$Whale))) { 
  df_i <- dplyr::filter(df, Whale == unique(df_2018$Whale)[i]) # We select one fin whale individual from 2018
  growth_rate <- 16  # Baleen plate growth rate (centimeters per year)
  last_sample <- unique(ymd(df_i$Data_capt)) # Date of last keratin sample (catch date)
  df_i <- 
    df_i %>%
    dplyr::mutate(days = Cm * 365.25 / growth_rate) %>% # Number of days between samples starting at Cm 0
    dplyr::mutate(days = days - days[1]) %>% # Number of days between samples starting at Cm -4
    dplyr::mutate(rev_days = days - days[length(days)]) %>% # Number of days starting at maximum Cm
    dplyr::mutate(sample_date = last_sample + rev_days) %>% # Estimated date
    dplyr::mutate(year = year(sample_date)) %>% # Estimated date year
    dplyr::mutate(year_rev = rev(sample_date)) # Estimated date retrospectively (real date)
  
  df_i$Year_from_sample_date <- format(as.Date(df_i$year_rev), "%Y") # year from the date estimated retrospectively
  df_i$Year_month <- format(as.Date(df_i$year_rev), "%Y/%m") # year and month from the date estimated retrospectively
  df_i$Month_from_sample_date <- format(as.Date(df_i$year_rev), "%m") # month from the date estimated retrospectively
  
  # plot top and bottom axis
  isop <- ggplot(df_i, aes(x = rev(Cm))) + 
    theme_classic(base_size = 16) + 
    xlab("Baleen sample (cm from youngest sample)") + 
    scale_x_continuous(breaks = df_i$Cm, 
                       labels = rev(df_i$Cm),
                       sec.axis = sec_axis(~., 
                                           breaks = rev(df_i$Cm),
                                           labels = rev(df_i$sample_date),
                                           name = "Date"))
  
  #isop <- ggplot(df_i, aes(x = rev(sample_date))) + xlab("Date") + scale_x_continuous(breaks = seq(min(df_i$sample_date), max(df_i$sample_date), by = "1 month"),labels = format(seq(min(rev(df_i$sample_date)), max(rev(df_i$sample_date)), by = "1 month"), "%b-%Y"))
  
  # plot carbon stable isotope values
  isop <- isop + geom_line(aes(y = dC), col = "black") + 
    geom_point(aes(y = dC), size = point_size, shape = 21, fill = "black", col = "black") +
    ylab(expression(paste(delta^{13}, "C (\u2030)")))+ theme_article()  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  a <- mean(na.omit(df_i$dN))# add the d15N data but need to scale it to d13C
  a
  b <- mean(na.omit(df_i$dC))# add the d15N data but need to scale it to d13C
  b
  v_shift <- a - b
  
  # plot nitrogen stable isotope values
  isop <- isop + geom_line(aes(y = dN - v_shift), col = "darkgrey") + 
    geom_point(aes(y = (dN - v_shift)), size = point_size, 
               shape = 21, fill = "grey", col = "darkgrey")
  
  # plot x and y axis
  isop <- isop + 
    scale_y_continuous(sec.axis = sec_axis(~.+v_shift, 
                                           name = expression(paste(delta^{15}, "N (\u2030)"))),limits = c(-30, -10)) + labs(title= unique(df_2018$Whale)[i])
  new_df <- df_i
  combined_df_2018 <- rbind(combined_df_2018, new_df)
  print(isop)
}

# plot time-calibrated stable isotope data for the individuals of 2018
nb.cols <- 13
colorRampPalette(brewer.pal(9, "Oranges"))(nb.cols) #palette for plot
unique(combined_df_2018$Whale)

plot_dN_2018 <- ggplot() + geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18044",], aes(x= rev(sample_date), y= dN),col="#FFF5EB")+ 
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18044",],aes(x= rev(sample_date), y= dN), size = point_size, shape = 21, fill = "#FFF5EB") +
  geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18043",], aes(x= rev(sample_date), y= dN),col="#FEEBD7")+ 
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18043",],aes(x= rev(sample_date), y= dN), size = point_size, shape = 21, fill = "#FEEBD7") + 
  geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18038",], aes(x= rev(sample_date), y= dN),col="#FDDEBF") +
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18038",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#FDDEBF") +
  geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18050",], aes(x= rev(sample_date), y= dN),col="#FDD0A2")+ 
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18050",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#FDD0A2") + 
  geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18051",], aes(x= rev(sample_date), y= dN),col="#FDB97D")+ 
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18051",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#FDB97D") +
  geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18030",], aes(x= rev(sample_date), y= dN),col="#FDA35B")+ 
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18030",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#FDA35B") +
  geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18054",], aes(x= rev(sample_date), y= dN),col="#FD8D3C")+ 
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18054",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#FD8D3C") +
  geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18007",], aes(x= rev(sample_date), y= dN),col="#F57520")+ 
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18007",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#F57520") +
  geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18003",], aes(x= rev(sample_date), y= dN),col="#E95E0D")+ 
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18003",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#E95E0D") +
  geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18009",], aes(x= rev(sample_date), y= dN),col="#D94801")+ 
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18009",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#D94801") +
  geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18025",], aes(x= rev(sample_date), y= dN),col="#B73C02")+ 
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18025",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#B73C02") +
  geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18098",], aes(x= rev(sample_date), y= dN),col="#993103")+ 
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18098",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#993103") +
  geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18022",], aes(x= rev(sample_date), y= dN),col="#7F2704")+ 
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18022",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#7F2704") +
  theme_article(base_size = 15) + 
  xlab("Year") + theme(aspect.ratio = 3/4) +
  ylab(expression(paste(delta^{15}, "N (\u2030)")))

plot_dC_2018 <- ggplot() + geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18044",], aes(x= rev(sample_date), y= dC),col="#FFF5EB")+ 
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18044",],aes(x= rev(sample_date), y= dC), size = point_size, shape = 21, fill = "#FFF5EB") +
  geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18043",], aes(x= rev(sample_date), y= dC),col="#FEEBD7")+ 
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18043",],aes(x= rev(sample_date), y= dC), size = point_size, shape = 21, fill = "#FEEBD7") + 
  geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18038",], aes(x= rev(sample_date), y= dC),col="#FDDEBF") +
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18038",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#FDDEBF") +
  geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18050",], aes(x= rev(sample_date), y= dC),col="#FDD0A2")+ 
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18050",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#FDD0A2") + 
  geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18051",], aes(x= rev(sample_date), y= dC),col="#FDB97D")+ 
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18051",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#FDB97D") +
  geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18030",], aes(x= rev(sample_date), y= dC),col="#FDA35B")+ 
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18030",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#FDA35B") +
  geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18054",], aes(x= rev(sample_date), y= dC),col="#FD8D3C")+ 
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18054",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#FD8D3C") +
  geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18007",], aes(x= rev(sample_date), y= dC),col="#F57520")+ 
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18007",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#F57520") +
  geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18003",], aes(x= rev(sample_date), y= dC),col="#E95E0D")+ 
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18003",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#E95E0D") +
  geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18009",], aes(x= rev(sample_date), y= dC),col="#D94801")+ 
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18009",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#D94801") +
  geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18025",], aes(x= rev(sample_date), y= dC),col="#B73C02")+ 
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18025",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#B73C02") +
  geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18098",], aes(x= rev(sample_date), y= dC),col="#993103")+ 
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18098",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#993103") +
  geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18022",], aes(x= rev(sample_date), y= dC),col="#7F2704")+ 
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18022",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#7F2704") +
  theme_article(base_size = 15) + 
  xlab("Year") + theme(aspect.ratio = 3/4) +
  ylab(expression(paste(delta^{13}, "C (\u2030)")))

plot_dN_2018/plot_dC_2018

# Average values for 2018 individuals

#Nitrogen stable isotope values

#merge_2018_filtered <- combined_df_2018 %>% dplyr::group_by(Cm) %>% dplyr::summarise(across(.cols = where(is.numeric), .fns = ~mean(dN, na.rm = TRUE))) # Average value per centimeter
#ggplot(merge_2018_filtered, aes(x=merge_2018_filtered$Cm,y=merge_2018_filtered$dN)) + geom_point() + geom_line() + theme_article() + xlab("Cm") + ylab(expression(paste(delta^{15}, "N (\u2030)")))  

#Carbon stable isotope values

#merge_2018_filtered <- combined_df_2018 %>% dplyr::group_by(Cm) %>% dplyr::summarise(across(.cols = where(is.numeric), .fns = ~mean(dC, na.rm = TRUE))) # Average value per centimeter
#ggplot(merge_2018_filtered, aes(x=merge_2018_filtered$Cm,y=merge_2018_filtered$dC)) + geom_point() + geom_line() + theme_article() + xlab("Cm") + ylab(expression(paste(delta^{13}, "C (\u2030)")))  

#________________________________________________________________________________
# Samples from 2022 whales
#________________________________________________________________________________

df_2022 <- dplyr::filter(df, Year == 2022)
unique(df_2022$Whale) # (n=4; "F22046" "F22053" "F22054" "F22056")

point_size <- 2
combined_df_2022 <- NULL;

for (i in 1:length(unique(df_2022$Whale))) { 
  df_i <- dplyr::filter(df, Whale == unique(df_2022$Whale)[i]) # We select one fin whale individual from 2022
  growth_rate <- 16  # Baleen plate growth rate (centimeters per year)
  last_sample <- unique(ymd(df_i$Data_capt)) # Date of last keratin sample (catch date)
  df_i <- 
    df_i %>%
    dplyr::mutate(days = Cm * 365.25 / growth_rate) %>% # Number of days between samples starting at Cm 0
    dplyr::mutate(days = days - days[1]) %>% # Number of days between samples starting at Cm -4
    dplyr::mutate(rev_days = days - days[length(days)]) %>% # Number of days starting at maximum Cm
    dplyr::mutate(sample_date = last_sample + rev_days) %>% # Estimated date
    dplyr::mutate(year = year(sample_date)) %>% # Estimated date year
    dplyr::mutate(year_rev = rev(sample_date)) # Estimated date retrospectively (real date)
  
  df_i$Year_from_sample_date <- format(as.Date(df_i$year_rev), "%Y") # year from the date estimated retrospectively
  df_i$Year_month <- format(as.Date(df_i$year_rev), "%Y/%m") # year and month from the date estimated retrospectively
  df_i$Month_from_sample_date <- format(as.Date(df_i$year_rev), "%m") # month from the date estimated retrospectively
  
  # plot top and bottom axis
  isop <- ggplot(df_i, aes(x = rev(Cm))) + 
    theme_classic(base_size = 16) + 
    xlab("Baleen sample (cm from youngest sample)") + 
    scale_x_continuous(breaks = df_i$Cm, 
                       labels = rev(df_i$Cm),
                       sec.axis = sec_axis(~., 
                                           breaks = rev(df_i$Cm),
                                           labels = rev(df_i$sample_date),
                                           name = "Date"))
  
  #isop <- ggplot(df_i, aes(x = rev(sample_date))) + xlab("Date") + scale_x_continuous(breaks = seq(min(df_i$sample_date), max(df_i$sample_date), by = "1 month"),labels = format(seq(min(rev(df_i$sample_date)), max(rev(df_i$sample_date)), by = "1 month"), "%b-%Y"))
  
  # plot carbon stable isotope values
  isop <- isop + geom_line(aes(y = dC), col = "black") + 
    geom_point(aes(y = dC), size = point_size, shape = 21, fill = "black", col = "black") +
    ylab(expression(paste(delta^{13}, "C (\u2030)")))+ theme_article()  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  a <- mean(na.omit(df_i$dN))# add the d15N data but need to scale it to d13C
  a
  b <- mean(na.omit(df_i$dC))# add the d15N data but need to scale it to d13C
  b
  v_shift <- a - b
  
  # plot nitrogen stable isotope values
  isop <- isop + geom_line(aes(y = dN - v_shift), col = "darkgrey") + 
    geom_point(aes(y = (dN - v_shift)), size = point_size, 
               shape = 21, fill = "grey", col = "darkgrey")
  
  # plot x and y axis
  isop <- isop + 
    scale_y_continuous(sec.axis = sec_axis(~.+v_shift, 
                                           name = expression(paste(delta^{15}, "N (\u2030)"))),limits = c(-30, -10)) + labs(title= unique(df_2022$Whale)[i])
  new_df <- df_i
  combined_df_2022 <- rbind(combined_df_2022, new_df)
  print(isop)
}

# plot time-calibrated stable isotope data for the individuals of 2022
brewer.pal(n = 4, name = "Blues") #palette for plot
unique(combined_df_2022$Whale)

plot_dN_2022 <- ggplot() + geom_line(data= combined_df_2022[combined_df_2022$Whale == "F22046",], aes(x= rev(sample_date), y= dN),col="#EFF3FF")+ 
  geom_point(data= combined_df_2022[combined_df_2022$Whale == "F22046",],aes(x= rev(sample_date), y= dN), size = point_size, shape = 21, fill = "#EFF3FF") +
  geom_line(data= combined_df_2022[combined_df_2022$Whale == "F22053",], aes(x= rev(sample_date), y= dN),col="#BDD7E7")+ 
  geom_point(data= combined_df_2022[combined_df_2022$Whale == "F22053",],aes(x= rev(sample_date), y= dN), size = point_size, shape = 21, fill = "#BDD7E7") + 
  geom_line(data= combined_df_2022[combined_df_2022$Whale == "F22054",], aes(x= rev(sample_date), y= dN),col="#6BAED6") +
  geom_point(data= combined_df_2022[combined_df_2022$Whale == "F22054",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#6BAED6") +
  geom_line(data= combined_df_2022[combined_df_2022$Whale == "F22056",], aes(x= rev(sample_date), y= dN),col="#2171B5")+ 
  geom_point(data= combined_df_2022[combined_df_2022$Whale == "F22056",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#2171B5") + 
  theme_article(base_size = 15) + 
  xlab("Year") + theme(aspect.ratio = 3/4) +
  ylab(expression(paste(delta^{15}, "N (\u2030)")))

plot_dC_2022 <- ggplot() + geom_line(data= combined_df_2022[combined_df_2022$Whale == "F22046",], aes(x= rev(sample_date), y= dC),col="#EFF3FF")+ 
  geom_point(data= combined_df_2022[combined_df_2022$Whale == "F22046",],aes(x= rev(sample_date), y= dC), size = point_size, shape = 21, fill = "#EFF3FF") +
  geom_line(data= combined_df_2022[combined_df_2022$Whale == "F22053",], aes(x= rev(sample_date), y= dC),col="#BDD7E7")+ 
  geom_point(data= combined_df_2022[combined_df_2022$Whale == "F22053",],aes(x= rev(sample_date), y= dC), size = point_size, shape = 21, fill = "#BDD7E7") + 
  geom_line(data= combined_df_2022[combined_df_2022$Whale == "F22054",], aes(x= rev(sample_date), y= dC),col="#6BAED6") +
  geom_point(data= combined_df_2022[combined_df_2022$Whale == "F22054",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#6BAED6") +
  geom_line(data= combined_df_2022[combined_df_2022$Whale == "F22056",], aes(x= rev(sample_date), y= dC),col="#2171B5")+ 
  geom_point(data= combined_df_2022[combined_df_2022$Whale == "F22056",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#2171B5") + 
  theme_article(base_size = 15) + 
  xlab("Year") + theme(aspect.ratio = 3/4) +
  ylab(expression(paste(delta^{13}, "C (\u2030)")))

plot_dN_2022/plot_dC_2022

# Average values for 2022 individuals

#Nitrogen stable isotope values

#merge_2022_filtered <- combined_df_2022 %>% dplyr::group_by(Cm) %>% dplyr::summarise(across(.cols = where(is.numeric), .fns = ~mean(dN, na.rm = TRUE))) # Average value per centimeter
#ggplot(merge_2022_filtered, aes(x=merge_2022_filtered$Cm,y=merge_2022_filtered$dN)) + geom_point() + geom_line() + theme_article() + xlab("Cm") + ylab(expression(paste(delta^{15}, "N (\u2030)")))  

#Carbon stable isotope values

#merge_2022_filtered <- combined_df_2022 %>% dplyr::group_by(Cm) %>% dplyr::summarise(across(.cols = where(is.numeric), .fns = ~mean(dC, na.rm = TRUE))) # Average value per centimeter
#ggplot(merge_2022_filtered, aes(x=merge_2022_filtered$Cm,y=merge_2022_filtered$dC)) + geom_point() + geom_line() + theme_article() + xlab("Cm") + ylab(expression(paste(delta^{13}, "C (\u2030)")))  

#__________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
# Merge the individuals from 2013,2015,2018 and 2022 data
#__________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

plot_iso_temp <- (plot_dN_2022+plot_dC_2022) / (plot_dN_2018+plot_dC_2018) / (plot_dN_2015+plot_dC_2015) / (plot_dN_2013+plot_dC_2013) 
ggsave("/Users/marcruizisagales/Documents/GitHub/Climate-baleen-plate-isotopes/png/Figure_S4_Isotopic_time_series_per_years.png", plot_iso_temp, 
       device = png(width = 400, height = 600))

merge <- rbind(combined_df_2013,combined_df_2015,combined_df_2018, combined_df_2022)

library(openxlsx) #Save
write.xlsx(merge, file = "/Users/marcruizisagales/Documents/GitHub/Climate-baleen-plate-isotopes/data/All_merged_time_calibrated_2013_to_2022_non_Suess_cor.xlsx")



