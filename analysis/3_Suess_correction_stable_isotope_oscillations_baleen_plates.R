#--------------------------------------------------------------------------------
#Suess correction of time-calibrated stable isotope data
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
library(dplyr)
library(SuessR)

# 2. Import data
merge <- read_excel("/Users/marcruizisagales/Documents/GitHub/Climate-baleen-plate-isotopes/data/All_merged_time_calibrated_2013_to_2022_non_Suess_cor.xlsx")

merge <-merge[complete.cases(merge$dC), ] # remove the rows with no d13C values in order to be able to correct for the Suess effect

merge1 <- merge %>% 
  dplyr::mutate(Species = ifelse(Whale %in% c("F13129", "F18022", "F18098"), "hybrid", "fin")) 

merge2 <- dplyr::filter(merge1, Species == "fin") # remove blue-fin hybrids

length(unique(merge2$Whale)) # we have 30 fin whale individuals (being non-hybrids)

merge3 <- dplyr::filter(merge2, Whale != "F18025") # select all the whales that are not F18025 (non well-sampled)



# 3. Suess and Laws effect correction

merge3

merge3$id <-seq.int(nrow(merge3)) # add a sequence number to each row

merge3 <- merge3 %>% add_column(region = "Subpolar North Atlantic") # add the Subpolar North Atlantic region to each row

df1 <- merge3 # rename merge3 to df1

subset <- df1[c("id", "dC","Year_from_sample_date","region")] # select the merge3 column names (for RSuess) and name it subset

names(subset) <- c("id", "d13c","year","region") # rename subset columnames

subset$year <- as.numeric(subset$year) # define year as a numeric variable

subset <- as.data.frame(subset) # define subset as a dataframe

df2 <- SuessR(data=subset, correct.to = 2022) # correct the Suess and Laws effect

range(df2$net.cor)

fw.data_d13cor <- merge(df1,df2,by="id") # merge df1 and df2

fw.data_d13cor_summary_whale <- fw.data_d13cor %>% dplyr::group_by(Whale) %>% dplyr::summarise(dN_mean= mean(dN, na.rm=T),
                                                 dN_sd= sd(dN, na.rm=T),
                                                 dN_min= min(dN, na.rm=T),
                                                 dN_max= max(dN, na.rm=T),
                                                 dN_max= max(dN, na.rm=T),
                                                 dN_range= max(dN, na.rm=T)-min(dN, na.rm=T),
                                                 dC_sd= sd(d13c.cor, na.rm=T),
                                                 dC_min= min(d13c.cor, na.rm=T),
                                                 dC_max= max(d13c.cor, na.rm=T),
                                                 dC_range= max(d13c.cor, na.rm=T)-min(d13c.cor, na.rm=T))

#overall
range(fw.data_d13cor$dN, na.rm=T)
range(fw.data_d13cor$d13c.cor)

#per whale
mean(fw.data_d13cor_summary_whale$dN_range)
sd(fw.data_d13cor_summary_whale$dN_range)

mean(fw.data_d13cor_summary_whale$dC_range)
sd(fw.data_d13cor_summary_whale$dC_range)

library(openxlsx) # Save
write.xlsx(fw.data_d13cor, file = "/Users/marcruizisagales/Documents/GitHub/Climate-baleen-plate-isotopes/data/All_merged_time_calibrated_2013_to_2022_Suess_cor.xlsx")
