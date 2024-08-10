#--------------------------------------------------------------------------------
#Sliding window analysis with the time-calibrated stable isotope data (nitrogen)
#--------------------------------------------------------------------------------

# adapted from DR de Zwaan, A Drake, JL Greenwood, and K Martin (2020).
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
library(plyr)
library(dplyr)
library(lme4)
library(climwin)
library(readr)

# 2. Import stable isotope data
df <- read_excel("/Users/marcruizisagales/Documents/GitHub/Climate-baleen-plate-isotopes/data/All_merged_time_calibrated_2013_to_2022_Suess_cor.xlsx") # Import Suess-corrected dataset of stable isotope data along the baleen plate of fin whales 

length(unique(df$Whale)) # check the 29 individuals
min(df$year_rev) # first date: "2010-07-11 UTC"
max(df$year_rev) # last date: "2022-07-30 UTC"

# 3. Load monthly climate data
data_climate <- read_excel("/Users/marcruizisagales/Desktop/Papers/Paper_climate/Datasets_paper_clim/merged_climate_cycles_dated_ts.xlsx") #load climate data
data_climate<- as.data.frame(data_climate)
data_climate$Date <- as.Date(data_climate$date)
data_climate$date_use <- format(data_climate$Date, format= "%d/%m/%Y")

data_env <- read_excel("/Users/marcruizisagales/Desktop/Papers/Paper_climate/Datasets_paper_clim/Env_total_13_nov.xlsx") # load environmental data
data_env<- as.data.frame(data_env)
data_env$Date <- as.Date(data_env$Date)

data_climate<- merge(data_climate, data_env, by = "Date", all=TRUE)

# 4. Convert measurement date into proper date format
df$measured <- strptime(paste(df$year_rev), format="%Y-%m-%d")
df$measured_date_formatted <- format(df$measured, format= "%d/%m/%Y")

# 5. Select only samples with isotopic measurements
Bp_isotopic_data_use <-df[complete.cases(df$dN), ]
Bp_isotopic_data_use <- Bp_isotopic_data_use[complete.cases(Bp_isotopic_data_use$d13c.cor), ]
nrow(Bp_isotopic_data_use)

range(Bp_isotopic_data_use$dN)
range(Bp_isotopic_data_use$d13c.cor)

mean(Bp_isotopic_data_use$dN)
mean(Bp_isotopic_data_use$d13c.cor)

sd(Bp_isotopic_data_use$dN)
sd(Bp_isotopic_data_use$d13c.cor)

# 6. Convert year and whale ID to factor
Bp_isotopic_data_use$month_f <- factor(Bp_isotopic_data_use$Month_from_sample_date)
Bp_isotopic_data_use$year_f <- factor(Bp_isotopic_data_use$year.y)
Bp_isotopic_data_use$whale_id <- factor(Bp_isotopic_data_use$Whale)

# 7. Make julian date (yearday)
Bp_isotopic_data_use$julian_day <-as.POSIXlt(Bp_isotopic_data_use$year_rev)$yday

# 8. Dataframe summary
length(unique(Bp_isotopic_data_use$whale_id)) # Number of whales = 29
length(unique(Bp_isotopic_data_use$year_f)) #Number of years = 13
Samples_per_whale <- plyr::count(Bp_isotopic_data_use, "whale_id") # Number of samples along the baleen plate per whale
Samples_per_whale
mean(Samples_per_whale$freq) # 47.65517 +/- 2.56732 samples per baleen
sd(Samples_per_whale$freq)
sum(Samples_per_whale$freq) # Total samples = 1381 

# 9. Add sample ID

Bp_isotopic_data_use$sample_ID <- seq.int(nrow(Bp_isotopic_data_use)) #1381

# 10. Build null model

# 10.1. Test if sex should be included to control for among sex differences

wo_sex_model <- lmer(dN ~ poly(julian_day,2) + (1| whale_id), data=Bp_isotopic_data_use)
sex_model <- lmer(dN ~ poly(julian_day,2) + as.factor(Bp_isotopic_data_use$Sex) + (1| whale_id), data=Bp_isotopic_data_use)

anova(wo_sex_model, sex_model, test="Chisq") # Decision: Evidence that sex does not improve the model. Do not include.

# 10.2. Test if lenght should be included to control for among lenght differences

wo_lenght_model <- lmer(dN ~ poly(julian_day,2) + (1| whale_id), data=Bp_isotopic_data_use)
lenght_model <- lmer(dN ~ poly(julian_day,2) + as.factor(Bp_isotopic_data_use$Length) + (1| whale_id), data=Bp_isotopic_data_use)

anova(wo_lenght_model, lenght_model, test="Chisq") # Decision: Evidence that lenght does not improve the model. Lenght almost significative.

# 10.3. Test if reproductive status should be included to control for among status differences

wo_status_model <- lmer(dN ~ poly(julian_day,2) + (1| whale_id), data=Bp_isotopic_data_use)
status_model <- lmer(dN ~ poly(julian_day,2) + as.factor(Bp_isotopic_data_use$Status) + (1| whale_id), data=Bp_isotopic_data_use)

anova(wo_status_model, status_model, test="Chisq") # Decision: Evidence that reproductive status does not improve the model. Do not include.

# 11. Sliding windows

Bp_dN_clim_m_rel <- 
  climwin::slidingwin(baseline = lme4::lmer(dN ~ 1 + poly(Bp_isotopic_data_use$julian_day,2) + (1| Bp_isotopic_data_use$whale_id), data=Bp_isotopic_data_use),  #this is the base model to which the program adds the climate variables
                      xvar = list(NAO=data_climate$NAO_index,
                                  AMO= data_climate$AMO_index,
                                  AMOC= data_climate$AMOC_index),
                      type = "relative", 
                      range= c(36,0), 
                      stat = c("mean"),
                      func = c("lin"),
                      cinterval = "month",
                      cmissing = "method2", # only applied in two months of AMOC index
                      cdate = data_climate$date_use, bdate = Bp_isotopic_data_use$measured_date_formatted)

Bp_dN_clim_m_rel$combos #check all of the sliding windows in the analyses

#AIC difference < 2
a<-head(Bp_dN_clim_m_rel[[1]]$Dataset, 20) ## Lin_func Mean NAO, close window-open window (34-31) (beta coefficient= -0.36)
a <- a %>% filter(deltaAICc < (a$deltaAICc[1])+2)
nrow_a <- nrow(a)
a$Index <- c(rep("NAO", nrow_a))
b<-head(Bp_dN_clim_m_rel[[2]]$Dataset, 20) ## Lin_func Mean AMO, close window-open window (32-27) (beta coefficient= 2.57)
b <- b %>% filter(deltaAICc < (b$deltaAICc[1])+2)
nrow_b <- nrow(b)
b$Index <- c(rep("AMO", nrow_b))
c<-head(Bp_dN_clim_m_rel[[3]]$Dataset, 20) ## Lin_func Mean AMOC, close window-open window (30-13) (beta coefficient= 0.35)
c <- c %>% filter(deltaAICc < (c$deltaAICc[1])+2)
nrow_c <- nrow(c)
c$Index <- c(rep("AMOC", nrow_c))

data <- rbind(a,b,c)

levels=c("NAO","AMO","AMOC")

data$Index <- factor(data$Index, levels=rev(c("NAO","AMO","AMOC")))

ggplot(data) +
  geom_vline(aes(xintercept = data$Index), color = "gray", linetype = "dotted") +
  geom_linerange(aes(x = data$Index, ymin = WindowOpen, ymax = WindowClose), color = '#B3CDE3', linewidth = 20, alpha = 0.5) +  # Afegeix les línies de separació
  geom_hline(yintercept = seq(1,36,1), lty = 1, lwd = 0.1, colour = "white") +
  geom_hline(yintercept = 0, lty = 1, lwd = 1, colour = "yellow") +
  coord_flip() +
  labs(x = "", y = "Time from last keratin sample (months)", title = "") +
  theme_article() + ylim(36,0)+
  scale_x_discrete(breaks=levels,labels=c("Avg NAO index","Avg AMO index","Avg AMOC index")) + theme(aspect.ratio=3/4)


# 12.Randomization
#Bp_dN_rand_clim_m_rel <- randwin(repeats = 1000, 
#baseline = lme4::lmer(dN ~ 1 + poly(Bp_isotopic_data_use$julian_day,2) + (1| Bp_isotopic_data_use$whale_id), data=Bp_isotopic_data_use),  #this is the base model to which the program adds the climate variables
#xvar = list(NAO=data_climate$NAO_index,
#AMO= data_climate$AMO_index,
#AMOC= data_climate$AMOC_index),
#type = "relative", 
#range= c(36,0), 
#stat = c("mean"),
#func = c("lin"),
#cinterval = "month",
#cmissing = "method2", # only applied in two months of AMOC index
#cdate = data_climate$date_use, bdate = Bp_isotopic_data_use$measured_date_formatted)

#Bp_dN_rand_clim_m_rel



# 13. Compare randomizations to original sliding window


# Mean NAO
pvalue(dataset = Bp_dN_clim_m_rel[[1]]$Dataset, datasetrand = Bp_dN_rand_clim_m_rel[[1]], 
       metric="AIC", sample.size=12) # "<0.001"

# Mean AMO
pvalue(dataset = Bp_dN_clim_m_rel[[2]]$Dataset, datasetrand = Bp_dN_rand_clim_m_rel[[2]], 
       metric="AIC", sample.size=12) # "<0.001"

# Mean AMOC
pvalue(dataset = Bp_dN_clim_m_rel[[3]]$Dataset, datasetrand = Bp_dN_rand_clim_m_rel[[3]], 
       metric="AIC", sample.size=12) # "<0.001"


# 14. Extract weather variables from time periods identified by sliding window

# NAO index
#Best window for mean NAO linear 34-31
Bp_isotopic_data_use$mean_NAO_linear_start_m <- ymd(Bp_isotopic_data_use$measured)%m-% months(34)
Bp_isotopic_data_use$mean_NAO_linear_end_m <- ymd(Bp_isotopic_data_use$measured)%m-% months(31)

# AMO index
#Best window for mean AO linear 32-27
Bp_isotopic_data_use$mean_AMO_linear_start_m <- ymd(Bp_isotopic_data_use$measured)%m-% months(32)
Bp_isotopic_data_use$mean_AMO_linear_end_m <- ymd(Bp_isotopic_data_use$measured)%m-% months(27)

# AMOC
#Best window for mean AMOC linear 30-14
Bp_isotopic_data_use$mean_AMOC_linear_start_m <- ymd(Bp_isotopic_data_use$measured)%m-% months(30)
Bp_isotopic_data_use$mean_AMOC_linear_end_m <- ymd(Bp_isotopic_data_use$measured)%m-% months(14)

# Select entire 36-month period for all climate variables for summary table
Bp_isotopic_data_use$overall_start <- ymd(Bp_isotopic_data_use$measured)%m-% months(36)
Bp_isotopic_data_use$overall_end <- ymd(Bp_isotopic_data_use$measured)%m-% months(0)

# Convert dates to proper POSIXct format

data_climate$date.int <- strptime(paste(data_climate$date), format="%Y-%m-%d") 
data_climate$date_new <- as.POSIXct(data_climate$date.int)
head(data_climate)
data_climate <- data_climate[,-18] #make sure we delete date_int column
data_climate


# 15. Format conversion
# Average NAO index date
Bp_isotopic_data_use$mean_NAO_start_m_int <- strptime(paste(Bp_isotopic_data_use$mean_NAO_linear_start_m), format="%Y-%m-%d") 
Bp_isotopic_data_use$mean_NAO_start_m_date <- as.POSIXct(Bp_isotopic_data_use$mean_NAO_start_m_int)

Bp_isotopic_data_use$mean_NAO_end_m_int <- strptime(paste(Bp_isotopic_data_use$mean_NAO_linear_end_m), format="%Y-%m-%d") 
Bp_isotopic_data_use$mean_NAO_end_m_date <- as.POSIXct(Bp_isotopic_data_use$mean_NAO_end_m_int)


# Average AMO index date
Bp_isotopic_data_use$mean_AMO_start_m_int <- strptime(paste(Bp_isotopic_data_use$mean_AMO_linear_start_m), format="%Y-%m-%d") 
Bp_isotopic_data_use$mean_AMO_start_m_date <- as.POSIXct(Bp_isotopic_data_use$mean_AMO_start_m_int)

Bp_isotopic_data_use$mean_AMO_end_m_int <- strptime(paste(Bp_isotopic_data_use$mean_AMO_linear_end_m), format="%Y-%m-%d") 
Bp_isotopic_data_use$mean_AMO_end_m_date <- as.POSIXct(Bp_isotopic_data_use$mean_AMO_end_m_int)


# Average AMOC index date
Bp_isotopic_data_use$mean_AMOC_start_m_int <- strptime(paste(Bp_isotopic_data_use$mean_AMOC_linear_start_m), format="%Y-%m-%d") 
Bp_isotopic_data_use$mean_AMOC_start_m_date <- as.POSIXct(Bp_isotopic_data_use$mean_AMOC_start_m_int)

Bp_isotopic_data_use$mean_AMOC_end_m_int <- strptime(paste(Bp_isotopic_data_use$mean_AMOC_linear_end_m), format="%Y-%m-%d") 
Bp_isotopic_data_use$mean_AMOC_end_m_date <- as.POSIXct(Bp_isotopic_data_use$mean_AMOC_end_m_int)

# Overall time period
Bp_isotopic_data_use$overall_start_int <- strptime(paste(Bp_isotopic_data_use$overall_start), format="%Y-%m-%d") 
Bp_isotopic_data_use$overall_start_date <- as.POSIXct(Bp_isotopic_data_use$overall_start_int)

Bp_isotopic_data_use$overall_end_int <- strptime(paste(Bp_isotopic_data_use$overall_end), format="%Y-%m-%d") 
Bp_isotopic_data_use$overall_end_date <- as.POSIXct(Bp_isotopic_data_use$overall_end_int)


# 16. Loop to extract climate indices for each sample

# Empty lists

mean_NAO_dN_list <- list()


mean_AMO_dN_list <- list()


mean_AMOC_dN_list <- list()


overall_list <- list()


for (i in 1:1381) { # loop from one to number of samples samples 
  # Average NAO index
  mean_NAO_dN_list[[i]] <- data_climate %>% 
    filter(date_new >= Bp_isotopic_data_use$mean_NAO_start_m_date[i] 
           & date_new <= Bp_isotopic_data_use$mean_NAO_end_m_date[i]) %>% 
    
    summarise(mean_NAO_dN = mean(NAO_index))
  
  
  # Average AMO index
  mean_AMO_dN_list[[i]] <- data_climate %>% 
    filter(date_new >= Bp_isotopic_data_use$mean_AMO_start_m_date[i] 
           & date_new <= Bp_isotopic_data_use$mean_AMO_end_m_date[i]) %>% 
    
    summarise(mean_AMO_dN = mean(AMO_index))
  
  
  # Average AMOC index
  mean_AMOC_dN_list[[i]] <- data_climate %>% 
    filter(date_new >= Bp_isotopic_data_use$mean_AMOC_start_m_date[i] 
           & date_new <= Bp_isotopic_data_use$mean_AMOC_end_m_date[i]) %>% 
    
    summarise(mean_AMOC_dN = mean(AMOC_index))
  
  # Overall
  overall_list[[i]] <- data_climate %>% 
    filter(date_new >= Bp_isotopic_data_use$overall_start_date[i] 
           & date_new <= Bp_isotopic_data_use$overall_end_date[i]) %>% 
    
    summarise(mean_NAO = mean(NAO_index, na.rm = TRUE),
              mean_AMO = mean(AMO_index, na.rm = TRUE),
              mean_AMOC = mean(AMOC_index, na.rm = TRUE))
  
  
  # Convert lists to rows per sample ID
  
  mean_NAO_dN_data <- bind_rows(mean_NAO_dN_list)
  
  mean_AMO_dN_data <- bind_rows(mean_AMO_dN_list)
  
  mean_AMOC_dN_data <- bind_rows(mean_AMOC_dN_list)
  
  overall_data <- bind_rows(overall_list)
  
  
  # Convert each variable to a column in a larger dataset
  Bp_dN_windows <- bind_cols(mean_NAO_dN_data,
                             mean_AMO_dN_data, 
                             mean_AMOC_dN_data,
                             overall_data)
  
}

head(Bp_dN_windows)


### Combine with nestquadg measurement data

head(Bp_isotopic_data_use)

Bp_isotopic_data_new <- Bp_isotopic_data_use # Extract only columns of interest from original data

Bp_dN_final <- dplyr::bind_cols(Bp_isotopic_data_new, Bp_dN_windows)


### Save data

write.csv(Bp_dN_final, "/Users/marcruizisagales/Documents/GitHub/Climate-baleen-plate-isotopes/data/Bp_isotopic_climate_summarised_final_nitrogen_mean_seasonality_for_Github.csv")

