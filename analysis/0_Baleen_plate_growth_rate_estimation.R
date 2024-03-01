#--------------------------------------------------------------------------------
# Baleen plate growth rate estimation
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
library(zoo)
library(TSA)

# 2. Import data

total <- read_excel("/Users/marcruizisagales/Documents/GitHub/Climate-baleen-plate-isotopes/data/0_point_alignments_29_juny_total_-4.xlsx") # Import dataset of the one-centimetre-spaced stable isotope data along the baleen plate of fin whales


# function to add text for the highest frequency spectral 
# densities to periodogram plots

periodogramText <- function(p, k){
  
  # extract frequencies and spectral densities from the object
  dd <- data.frame(freq = p$freq, spec = p$spec)
  
  # sort them into rank order
  ord <- dd[order(-dd$spec), ]
  
  # extract the top two densities
  top2  <-  head(ord, 15)
  
  # convert to distance along our sampling transect
  distance <- 1/top2$freq
  
  # add either the first or second
  text(0.4, max(p$spec) * 0.8, 
       paste0(sprintf("%.1f", round(distance[k], 1))," cm"), 
       cex = 2)
}

#Function to find max lag in ccf

Find_Max_CCF<- function(a,b)
{
  d <- ccf(a, b, plot = FALSE)
  cor = d$acf[,,1]
  lag = d$lag[,,1]
  res = data.frame(cor,lag)
  res_max = res[which.max(res$cor),]
  return(res_max)
} 


#2013
# Whale F13065= = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
whaleF13065 <- filter(total, Whale == "F13065")
par(mfrow = c(1, 1))
p1 <- periodogram(na.approx(whaleF13065$dN)) #k=1
#p1 <- periodogram(na.approx(whaleF13065$dC)) #k=1
periodogramText(p1, k = 1) #periode de 15 cm

ccf_p1 <- ccf(whaleF13065$dC, whaleF13065$dN, lag.max = 10, plot = TRUE)
ccf_p1_text <- Find_Max_CCF(whaleF13065$dC, whaleF13065$dN) #cor 0.721362 and lag 0
ccf_p1_text 
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Whale F13066= = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
whaleF13066 <- filter(total, Whale == "F13066")
par(mfrow = c(1, 1))
p2 <- periodogram(na.approx(whaleF13066$dN)) #k=1
#p2 <- periodogram(na.approx(whaleF13066$dC)) #k=1
periodogramText(p2, k = 1) #periode de 15 cm

ccf_p2 <- ccf(whaleF13066$dC, whaleF13066$dN, lag.max = 10, plot = TRUE)
ccf_p2_text <- Find_Max_CCF(whaleF13066$dC, whaleF13066$dN) #cor 0.661634 and lag 3
ccf_p2_text 
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Whale F13068= = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
whaleF13068 <- filter(total, Whale == "F13068")
par(mfrow = c(1, 1))
p3 <- periodogram(na.approx(whaleF13068$dN)) #k=1
#p3 <- periodogram(na.approx(whaleF13068$dC)) #k=4
periodogramText(p3, k = 1) #periode de 15 cm

ccf_p3 <- ccf(whaleF13068$dC, whaleF13068$dN, lag.max = 10, plot = TRUE)
ccf_p3_text <- Find_Max_CCF(whaleF13068$dC, whaleF13068$dN) #cor 0.5142207 and lag 7
ccf_p3_text
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Whale F13073= = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
whaleF13073 <- filter(total, Whale == "F13073")
par(mfrow = c(1, 1))
p4 <- periodogram(na.approx(whaleF13073$dN)) #k=1
#p4 <- periodogram(na.approx(whaleF13073$dC)) #k=1
p4_text <- periodogramText(p4, k = 1) #periode de 15 cm

ccf_p4 <- ccf(whaleF13073$dC, whaleF13073$dN, lag.max = 10, plot = TRUE)
ccf_p4_text <- Find_Max_CCF(whaleF13073$dC, whaleF13073$dN) #cor 0.5062958 and lag 3
ccf_p4_text
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Whale F13076= = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
whaleF13076 <- filter(total, Whale == "F13076")
par(mfrow = c(1, 1))
p5 <- periodogram(na.approx(whaleF13076$dN)) #k=1
#p5 <- periodogram(na.approx(whaleF13076$dC)) #k=2
p5_text <- periodogramText(p5, k = 1) #periode de 15 cm

ccf_p5 <- ccf(whaleF13076$dC, whaleF13076$dN, lag.max = 10, plot = TRUE)
ccf_p5_text <- Find_Max_CCF(whaleF13076$dC, whaleF13076$dN) #cor 0.5759791 and lag -3
ccf_p5_text
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

#2015
# Whale F15078= = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
whaleF15078 <- filter(total, Whale == "F15078")
par(mfrow = c(1, 1))
p6 <- periodogram(na.approx(whaleF15078$dN))  #k=2
#p6 <- periodogram(na.approx(whaleF15078$dC)) #k=1
periodogramText(p6, k = 2) #periode de 15 cm

ccf_p6 <- ccf(whaleF15078$dC, whaleF15078$dN, lag.max = 10, plot = TRUE)
ccf_p6_text <- Find_Max_CCF(whaleF15078$dC, whaleF15078$dN) #cor 0.6612612 and lag 0
ccf_p6_text 
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Whale F15097= = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
whaleF15097 <- filter(total, Whale == "F15097")
par(mfrow = c(1, 1))
p7 <- periodogram(na.approx(whaleF15097$dN)) #k=1
#p7 <- periodogram(na.approx(whaleF15097$dC)) #k=1
periodogramText(p7, k = 1) #periode de 16 cm

ccf_p7 <- ccf(whaleF15097$dC, whaleF15097$dN, lag.max = 10, plot = TRUE)
ccf_p7_text <- Find_Max_CCF(whaleF15097$dC, whaleF15097$dN) #cor 0.7150837 and lag 0
ccf_p7_text 
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Whale F15086= = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
whaleF15086 <- filter(total, Whale == "F15086")
par(mfrow = c(1, 1))
p8 <- periodogram(na.approx(whaleF15086$dN)) #k=2
#p8 <- periodogram(na.approx(whaleF15086$dC)) #k=3
periodogramText(p8, k = 2) #periode de 16.7 cm

ccf_p8 <- ccf(whaleF15086$dC, whaleF15086$dN, lag.max = 10, plot = TRUE)
ccf_p8_text <- Find_Max_CCF(whaleF15086$dC, whaleF15086$dN) #cor 0.554663 and lag -11
ccf_p8_text 
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Whale F15083= = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
whaleF15083 <- filter(total, Whale == "F15083")
par(mfrow = c(1, 1))
p9 <- periodogram(na.approx(whaleF15083$dN)) #k=2
#p9 <- periodogram(na.approx(whaleF15083$dC)) #k=1
periodogramText(p9, k = 2) #periode de 15 cm

ccf_p9 <- ccf(whaleF15083$dC, whaleF15083$dN, lag.max = 10, plot = TRUE)
ccf_p9_text <- Find_Max_CCF(whaleF15083$dC, whaleF15083$dN) #cor 0.672713 and lag 0
ccf_p9_text 
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Whale F15080= = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
whaleF15080 <- filter(total, Whale == "F15080")
par(mfrow = c(1, 1))
#p10 <- periodogram(na.approx(whaleF15080$dN)) #k=5
p10 <- periodogram(na.approx(whaleF15080$dC)) #k=1
#periodogramText(p10, k = 5) #periode de 16 cm amb nitrogen k=4
periodogramText(p10, k = 1) #periode de 16 cm amb carboni k=1

ccf_p10 <- ccf(whaleF15080$dC, whaleF15080$dN, lag.max = 10, plot = TRUE)
ccf_p10_text <- Find_Max_CCF(whaleF15080$dC, whaleF15080$dN) #cor 0.461124 and lag 2
ccf_p10_text 
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Whale F15079 (Lactating) = = = = = = = = = = = = = = = = = = = = = = = = 
whaleF15079 <- filter(total, Whale == "F15079")
par(mfrow = c(1, 1))
#p11 <- periodogram(na.approx(whaleF15079$dN)) # k=1
p11 <- periodogram(na.approx(whaleF15079$dC)) #k=5
#periodogramText(p11, k = 1) #periode de 18 cm amn nit k=1
periodogramText(p11, k = 5) #periode de 18 cm amn carb k=5

ccf_p11 <- ccf(whaleF15079$dC, whaleF15079$dN, lag.max = 10, plot = TRUE)
ccf_p11_text <- Find_Max_CCF(whaleF15079$dC, whaleF15079$dN) #cor 0.324277 and lag 5
ccf_p11_text 
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Whale F15084= = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
whaleF15084 <- filter(total, Whale == "F15084")
par(mfrow = c(1, 1))
#p12 <- periodogram(na.approx(whaleF15084$dN)) #k=4
p12 <- periodogram(na.approx(whaleF15084$dC)) #k=5
#periodogramText(p12, k = 4) #periode de 16.7 cm en nitrogen k=4
periodogramText(p12, k = 5) #periode de 16.7 cm en nitrogen k=5

ccf_p12 <- ccf(na.approx(whaleF15084$dC), na.approx(whaleF15084$dN), lag.max = 10, plot = TRUE)
ccf_p12_text <- Find_Max_CCF(na.approx(whaleF15084$dC), na.approx(whaleF15084$dN)) #cor 0.6025646 and lag 1
ccf_p12_text 
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Whale F15087= = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
whaleF15087 <- filter(total, Whale == "F15087")
par(mfrow = c(1, 1))
p13 <- periodogram(na.approx(whaleF15087$dN)) #k=2
#p13 <- periodogram(na.approx(whaleF15087$dC)) #k=10
periodogramText(p13, k = 2) #periode de 16 cm nit k=2
periodogramText(p13, k = 10) #periode de 16 cm nit k=10

ccf_p13 <- ccf(na.approx(whaleF15087$dC), na.approx(whaleF15087$dN), lag.max = 10, plot = TRUE)
ccf_p13_text <- Find_Max_CCF(na.approx(whaleF15087$dC), na.approx(whaleF15087$dN)) #cor 0.2763756 and lag -7
ccf_p13_text 
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Whale F15088= = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
whaleF15088 <- filter(total, Whale == "F15088")
par(mfrow = c(1, 1))
p14 <- periodogram(na.approx(whaleF15088$dN)) #k=3
#p14 <- periodogram(na.approx(whaleF15088$dC)) 
periodogramText(p14, k = 3) #periode de 16 cm K=3 nit
periodogramText(p14, k = 2) #periode de 16 cm k=2 car

ccf_p14 <- ccf(whaleF15088$dC, whaleF15088$dN, lag.max = 10, plot = TRUE)
ccf_p14_text <- Find_Max_CCF(whaleF15088$dC, whaleF15088$dN) #cor 0.4855914 and lag 0
ccf_p14_text 
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

#2018
# Whale F18044= = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
whaleF18044 <- filter(total, Whale == "F18044")
par(mfrow = c(1, 1))
p15 <- periodogram(na.approx(whaleF18044$dN)) #k=1
#p15 <- periodogram(na.approx(whaleF18044$dC)) #k=3
periodogramText(p15, k = 1) #periode de 16.7 cm k=1 niy
periodogramText(p15, k = 3) #periode de 16.7 cm k=3 carb

ccf_p15 <- ccf(whaleF18044$dC, whaleF18044$dN, lag.max = 10, plot = TRUE)
ccf_p15_text <- Find_Max_CCF(whaleF18044$dC, whaleF18044$dN) #cor 0.5017228 and lag 5
ccf_p15_text 
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Whale F18043= = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
whaleF18043 <- filter(total, Whale == "F18043")
par(mfrow = c(1, 1))
p16 <- periodogram(na.approx(whaleF18043$dN)) #k=7
#p16 <- periodogram(na.approx(whaleF18043$dC)) #k=1
periodogramText(p16, k = 7) #periode de 16.7 cm

ccf_p16 <- ccf(whaleF18043$dC, whaleF18043$dN, lag.max = 10, plot = TRUE)
ccf_p16_text <- Find_Max_CCF(whaleF18043$dC, whaleF18043$dN) #cor 0.4639773 and lag 1
ccf_p16_text 
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Whale F18038= = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
whaleF18038 <- filter(total, Whale == "F18038")
par(mfrow = c(1, 1))
p17 <- periodogram(na.approx(whaleF18038$dN)) #k=2
#p17 <- periodogram(na.approx(whaleF18038$dC)) #k=4
periodogramText(p17, k = 2) #periode de 16.7 cm

ccf_p17 <- ccf(whaleF18038$dC, whaleF18038$dN, lag.max = 10, plot = TRUE)
ccf_p17_text <- Find_Max_CCF(whaleF18038$dC, whaleF18038$dN) #cor 0.5639275 and lag -11
ccf_p17_text 
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Whale F18050= = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
whaleF18050 <- filter(total, Whale == "F18050")
par(mfrow = c(1, 1))
#p18 <- periodogram(na.approx(whaleF18050$dN)) #en aquest cas el nitrogen no funciona bé
p18 <- periodogram(na.approx(whaleF18050$dC)) #k=1
periodogramText(p18, k = 1) #periode de 16.7 cm

ccf_p18 <- ccf(whaleF18050$dC, whaleF18050$dN, lag.max = 10, plot = TRUE)
ccf_p18_text <- Find_Max_CCF(whaleF18050$dC, whaleF18050$dN) #cor 0.5639275 and lag -11
ccf_p18_text 
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Whale F18051= = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
whaleF18051 <- filter(total, Whale == "F18051")
par(mfrow = c(1, 1))
p19 <- periodogram(na.approx(whaleF18051$dN)) #k=1
#p19 <- periodogram(na.approx(whaleF18051$dC)) #k=3
periodogramText(p19, k = 1) #periode de 16.7 cm

ccf_p19 <- ccf(whaleF18051$dC, whaleF18051$dN, lag.max = 10, plot = TRUE)
ccf_p19_text <- Find_Max_CCF(whaleF18051$dC, whaleF18051$dN) #cor 0.3157882 and lag -8
ccf_p19_text 
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Whale F18030= = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
whaleF18030 <- filter(total, Whale == "F18030")
par(mfrow = c(1, 1))
p20 <- periodogram(na.approx(whaleF18030$dN)) #k=1
#p20 <- periodogram(na.approx(whaleF18030$dC)) #k=1
periodogramText(p20, k = 1) #periode de 16.7 cm

ccf_p20 <- ccf(whaleF18030$dC, whaleF18030$dN, lag.max = 10, plot = TRUE)
ccf_p20_text <- Find_Max_CCF(whaleF18030$dC, whaleF18030$dN) #cor 0.3872712 and lag -12
ccf_p20_text 
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Whale F18054= = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
whaleF18054 <- filter(total, Whale == "F18054")
par(mfrow = c(1, 1))
p21 <- periodogram(na.approx(whaleF18054$dN)) #k=3
#p21 <- periodogram(na.approx(whaleF18054$dC)) #k=5
periodogramText(p21, k = 3) #periode de 16.7 cm

ccf_p21 <- ccf(whaleF18054$dC, whaleF18054$dN, lag.max = 10, plot = TRUE)
ccf_p21_text <- Find_Max_CCF(whaleF18054$dC, whaleF18054$dN) #cor 0.5185896 and lag 1
ccf_p21_text 
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Whale F18007= = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
whaleF18007 <- filter(total, Whale == "F18007")
par(mfrow = c(1, 1))
p22 <- periodogram(na.approx(whaleF18007$dN)) #k=2
#p22 <- periodogram(na.approx(whaleF18007$dC)) #k=3
periodogramText(p22, k = 2) #periode de 16 cm

ccf_p22 <- ccf(na.approx(whaleF18007$dC), na.approx(whaleF18007$dN), lag.max = 10, plot = TRUE)
ccf_p22_text <- Find_Max_CCF(na.approx(whaleF18007$dC), na.approx(whaleF18007$dN)) #cor 0.3243816 and lag 1
ccf_p22_text 
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Whale F18003= = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
whaleF18003 <- filter(total, Whale == "F18003")
par(mfrow = c(1, 1))
p23 <- periodogram(na.approx(whaleF18003$dN)) #k=1 12.5 50 10 4.2 25 3.8 5 8.3 5.6 3.1
#p23 <- periodogram(na.approx(whaleF18003$dC)) #k=1 12.5 25 5 50 6.2 8.3 5.6 4.5 16.7 2
periodogramText(p23, k = 1) #periode de 12.5 cm

ccf_p23 <- ccf(na.approx(whaleF18003$dC), na.approx(whaleF18003$dN), lag.max = 10, plot = TRUE)
ccf_p23_text <- Find_Max_CCF(na.approx(whaleF18003$dC), na.approx(whaleF18003$dN)) #cor 0.6496165 and lag 0
ccf_p23_text 
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Whale F18009= = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
whaleF18009 <- filter(total, Whale == "F18009")
par(mfrow = c(1, 1))
p24 <- periodogram(na.approx(whaleF18009$dN)) #k=1
#p24 <- periodogram(na.approx(whaleF18009$dC)) #k=1
periodogramText(p24, k = 1) #periode de 16 cm

ccf_p24 <- ccf(na.approx(whaleF18009$dC), na.approx(whaleF18009$dN), lag.max = 10, plot = TRUE)
ccf_p24_text <- Find_Max_CCF(na.approx(whaleF18009$dC), na.approx(whaleF18009$dN)) #cor 0.5746649 and lag 0
ccf_p24_text 
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Whale F18025= = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
whaleF18025 <- filter(total, Whale == "F18025")
par(mfrow = c(1, 1))
p25 <- periodogram(na.approx(whaleF18025$dN)) #k=1
p25 <- periodogram(na.approx(whaleF18025$dC)) 
periodogramText(p25, k = 1) #periode de 24 cm

ccf_p25 <- ccf(na.approx(whaleF18025$dC), na.approx(whaleF18025$dN), lag.max = 10, plot = TRUE)
ccf_p25_text <- Find_Max_CCF(na.approx(whaleF18025$dC), na.approx(whaleF18025$dN)) #cor 0.3907993 and lag 1
ccf_p25_text 
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Whale F18098 (hybrid)= = = = = = = = = = = = = = = = = = = = = = = = = = 
whaleF18098 <- filter(total, Whale == "F18098")
par(mfrow = c(1, 1))
p26 <- periodogram(na.approx(whaleF18098$dN)) #k=1
#p26 <- periodogram(na.approx(whaleF18098$dC)) #k=2
periodogramText(p26, k = 1) #periode de 18 cm

ccf_p26 <- ccf(na.approx(whaleF18098$dC), na.approx(whaleF18098$dN), lag.max = 10, plot = TRUE)
ccf_p26_text <- Find_Max_CCF(na.approx(whaleF18098$dC), na.approx(whaleF18098$dN)) #cor 0.2284071 and lag 0
ccf_p26_text 
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Whale F18022 (hybrid)= = = = = = = = = = = = = = = = = = = = = = = = = = 
whaleF18022 <- filter(total, Whale == "F18022")
par(mfrow = c(1, 1))
p27 <- periodogram(na.approx(whaleF18022$dN)) #k=1
#p27 <- periodogram(na.approx(whaleF18022$dC)) # k=7
periodogramText(p27, k = 1) #periode de 18 cm

ccf_p27 <- ccf(na.approx(whaleF18022$dC), na.approx(whaleF18022$dN), lag.max = 10, plot = TRUE)
ccf_p27_text <- Find_Max_CCF(na.approx(whaleF18022$dC), na.approx(whaleF18022$dN)) #cor 0.4109479 and lag 0
ccf_p27_text 
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

#2022
# Whale F22046= = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
whaleF22046 <- filter(total, Whale == "F22046")
par(mfrow = c(1, 1))
p28 <- periodogram(na.approx(whaleF22046$dN)) #k=1
#p28 <- periodogram(na.approx(whaleF22046$dC)) #k=1
periodogramText(p28, k = 1) #periode de 16.7 cm

ccf_p28 <- ccf(na.approx(whaleF22046$dC), na.approx(whaleF22046$dN), lag.max = 10, plot = TRUE)
ccf_p28_text <- Find_Max_CCF(na.approx(whaleF22046$dC), na.approx(whaleF22046$dN)) #cor 0.6291124 and lag 1
ccf_p28_text 
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Whale F22053= = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
whaleF22053 <- filter(total, Whale == "F22053")
par(mfrow = c(1, 1))
p29 <- periodogram(na.approx(whaleF22053$dN)) #k=1
#p29 <- periodogram(na.approx(whaleF22053$dC)) #k=1
periodogramText(p29, k = 1) #periode de 16.7 cm

ccf_p29 <- ccf(na.approx(whaleF22053$dC), na.approx(whaleF22053$dN), lag.max = 10, plot = TRUE)
ccf_p29_text <- Find_Max_CCF(na.approx(whaleF22053$dC), na.approx(whaleF22053$dN)) #cor 0.6066229 and lag 0
ccf_p29_text 
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Whale F22054= = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
whaleF22054 <- filter(total, Whale == "F22054")
par(mfrow = c(1, 1))
p30 <- periodogram(na.approx(whaleF22054$dN)) #k=1
#p30 <- periodogram(na.approx(whaleF22054$dC)) #k=1
periodogramText(p30, k = 1) #periode de 16.7 cm

ccf_p30 <- ccf(na.approx(whaleF22054$dC), na.approx(whaleF22054$dN), lag.max = 10, plot = TRUE)
ccf_p30_text <- Find_Max_CCF(na.approx(whaleF22054$dC), na.approx(whaleF22054$dN)) #cor 0.5866235 and lag 2
ccf_p30_text 
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Whale F22056 (immadur)= = = = = = = = = = = = = = = = = = = = = = = = = 
whaleF22056 <- filter(total, Whale == "F22056")
par(mfrow = c(1, 1))
#p31 <- periodogram(na.approx(whaleF22056$dN)) 
p31 <- periodogram(na.approx(whaleF22056$dC)) 
periodogramText(p31, k = 1) # 54 6 27 9 10.8 3.9 7.7 13.5 4.5 4.9
periodogramText(p31, k = 1) # 54 10.8 18 7.7 3.4 6 2.2 5.4 2.7 9

ccf_p31 <- ccf(na.approx(whaleF22056$dC), na.approx(whaleF22056$dN), lag.max = 10, plot = TRUE)
ccf_p31_text <- Find_Max_CCF(na.approx(whaleF22056$dC), na.approx(whaleF22056$dN)) #cor 0.5490016 and lag 0
ccf_p31_text 
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


#mean growth rate = 16.07097 ± 2.459565 cm/year
mean(c(15,15,15,15,15,
       15,16,16.7,15,16,18,16.7,16,16,
       16.7, 16.7,16.7, 16.7,16.7, 16.7, 16.7, 16, 12.5, 16,24, 18,18,
       16.7, 16,7,16.7)) #=16.07097 (last immature individual not counted)
sd(c(15,15,15,15,15,
     15,16,16.7,15,16,18,16.7,16,16,
     16.7, 16.7,16.7, 16.7,16.7, 16.7, 16.7, 16, 12.5, 16,24, 18,18,
     16.7, 16,7,16.7)) #=2.459565
