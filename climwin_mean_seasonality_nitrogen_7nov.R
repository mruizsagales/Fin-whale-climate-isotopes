#Statistical Analysis - Paper 1 (Baleen plates and climate Marc Ruiz-Sagalés) WITH MEAN

#8 novembre 2023

##########################################################################################
##########################################################################################
### (1) Sliding window analysis R code from:

### "Timing and intensity of weather events shape nestling development strategies 
###  in three alpine breeding songbirds"

###  DR de Zwaan, A Drake, JL Greenwood, and K Martin (2020)
###  Frontiers in Ecology and Evolution	

### doi: 10.3389/fevo.2020.570034

### Run with R version 3.6.3
########################################################################################## 

### Set working directory

#setwd("C:/PATH")

### Load required packages

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

setwd("~/Desktop/Doctorat/Analisis_isotops_barbes/Projecte_barbes_clima/Environmental_vars/")
### Load isotopic baleen plate data

fw.data_d13cor <- read_excel("/Users/marcruizisagales/All_merged_2013_to_2022_Suess_20_oct.xlsx")
#unique(fw.data_d13cor$Whale) #check individuals

max(fw.data_d13cor$year_rev) #les dades arriben fins al "2022-07-30 UTC"

### Load weather data

#a) daily (cobreix la serie temporal isotopica)

#data_ao <- read_csv("Desktop/norm.daily.ao.cdas.z1000.19500101_current.csv")
#data_ao$date <- as.Date(with(data_ao, paste(year, month, day,sep="/")), "%Y/%m/%d")

#data_nao <- read_csv("Desktop/norm.daily.nao.cdas.z500.19500101_current.csv")
#data_nao$date <- as.Date(with(data_nao, paste(year, month, day,sep="/")), "%Y/%m/%d")

#put all data frames into list
#df_list_d <- list(data_ao, data_nao)      

#merge all data frames together
#env_data_d <- df_list_d %>% reduce(full_join, by='date')

#b) monthly 

## Load weather data (cobreix la serie temporal isotopica)

data_weather <- read_excel("~/Desktop/merged_climate_cycles_dated_ts.xlsx")
data_weather<- as.data.frame(data_weather)
data_weather$Date <- as.Date(data_weather$date)

env_data <- read_excel("~/Desktop/Env_total_13_nov.xlsx")
env_data<- as.data.frame(env_data)
env_data$Date <- as.Date(env_data$Date)

#env_data$Month <- month(env_data$Date)
#env_data$Year <- year(env_data$Date)
#env_data$ym_date <-as.Date(paste(as.numeric(env_data$Year), as.numeric(env_data$Month), "01", sep="-"))
#names(env_data)
#env_data <- env_data %>% dplyr::group_by(ym_date) %>% dplyr::summarise(thetao_atl = mean(thetao_atl),chl_ice = mean(chl_ice),thetao_atl = mean(thetao_atl),chl_atl = mean(chl_atl))

#env_data$Date <- env_data$ym_date
#env_data <- as.data.frame(env_data)

#comprovar que trec les columnes correctes
###########################################################################################
### (1) Sp: Balaenoptera physalus
###########################################################################################

#### Convert measurement date into proper date format

fw.data_d13cor$measured <- strptime(paste(fw.data_d13cor$year_rev), format="%Y-%m-%d")
#fw.data_d13cor <- fw.data_d13cor %>% filter(measured < '2021-06-30') #select only samples until this date
fw.data_d13cor$measured_date_formatted <- format(fw.data_d13cor$measured, format= "%d/%m/%Y")


### Select only samples with isotopic measurements

Bp_isotopic_data_use <-fw.data_d13cor[complete.cases(fw.data_d13cor$dN), ]
Bp_isotopic_data_use <-Bp_isotopic_data_use[complete.cases(Bp_isotopic_data_use$d13c.cor), ]

### Convert year and whale ID to factor
Bp_isotopic_data_use$month_f <- factor(Bp_isotopic_data_use$Month_from_sample_date)
Bp_isotopic_data_use$year_f <- factor(Bp_isotopic_data_use$year.y)
Bp_isotopic_data_use$whale_id <- factor(Bp_isotopic_data_use$Whale)

### Make julian date
Bp_isotopic_data_use$julian_day <-as.POSIXlt(Bp_isotopic_data_use$year_rev)$yday

### Bp sample size

length(unique(Bp_isotopic_data_use$whale_id)) # Number of whales = 29
length(unique(Bp_isotopic_data_use$year_f)) #Number of years = 13
Samples_per_whale <- plyr::count(Bp_isotopic_data_use, "whale_id") # Number of samples along the baleen plate per whale
Samples_per_whale
mean(Samples_per_whale$freq)
sd(Samples_per_whale$freq)
sum(Samples_per_whale$freq) # Total samples = 1383 (some have been removed)


##################
### Incorporate weather data 
##################


#tot junt
all_env_wea_data<- merge(env_data, data_weather, by = "Date", all=TRUE)

#save in excel file the chl values with atlantic extent 
library(writexl)
writexl::write_xlsx(all_env_wea_data, '/Users/marcruizisagales/Desktop/all_env_wea_data.xlsx')

all_env_wea_data$date_use <- format(all_env_wea_data$Date, format= "%d/%m/%Y")


### Add sample ID

Bp_isotopic_data_use$sample_ID <- seq.int(nrow(Bp_isotopic_data_use)) #1328

############################################################################################
############################# Balaenoptera physalus sliding window ###################################
############################################################################################

### Background example:
### Sliding window set to maximum of 30 days prior to measurement date (7-days post-hatch).
### For an average nest (4 eggs, 1 egg laid a day, incubation on penultimate egg, and
### 12 days incubation), this means earliest time is 9 days prior to the first egg laid.
### While it is possible weather events prior to arrival at the breeding site may influence
### maternal condition and offspring development, addressing this is not possible with the 
### available data. The point of this analysis is to test the influence of local weather 
### events on the breeding grounds that influence breeding decisions and offspring 
### development. Additionally, egg formation (when nutrients and hormones are deposited
### in the egg) occurs about 4-5 days prior to the lay date, so we believe 9 days should
### capture this period adequately. 


### Standardized steps for identifying important time periods:
# For each weather variable and trait or isotope, the top period >= 3 days is selected by AIC. #min and window max
# If several top periods differ by less than 2 AIC, then the period with the strongest Beta 
# coefficient is selected.
# If there are radically different time periods within the top models, each is selected
# and competed against each other in Step 2: Model fitting and selection.

### See linked paper for more details:
### doi: 10.3389/fevo.2020.570034


###############
### (A) dN
###############


### Test if year should be included to control for among year differences (REMOVE YEAR)

wo_year_model <- lmer(dN ~ poly(julian_day,2) + (1| whale_id), data=Bp_isotopic_data_use)
year_model <- lmer(dN ~ poly(julian_day,2) + year_f + (1| whale_id), data=Bp_isotopic_data_use)

anova(wo_year_model, year_model, test="Chisq")

### Decision: Evidence that year improves the model. (Pero no incloem pq volem veure diferencies entre anys).

### Test if sex should be included to control for among sex differences

wo_sex_model <- lmer(dN ~ poly(julian_day,2) + (1| whale_id), data=Bp_isotopic_data_use)
sex_model <- lmer(dN ~ poly(julian_day,2) + as.factor(Bp_isotopic_data_use$Sex) + (1| whale_id), data=Bp_isotopic_data_use)

anova(wo_sex_model, sex_model, test="Chisq")

### Decision: Evidence that sex does not improve the model. Do not include.

### Test if lenght should be included to control for among lenght differences

wo_lenght_model <- lmer(dN ~ poly(julian_day,2) + (1| whale_id), data=Bp_isotopic_data_use)
lenght_model <- lmer(dN ~ poly(julian_day,2) + as.factor(Bp_isotopic_data_use$Length) + (1| whale_id), data=Bp_isotopic_data_use)

anova(wo_lenght_model, lenght_model, test="Chisq")

### Decision: Evidence that lenght does not improve the model. Lenght casi significativa.

### Test if status should be included to control for among status differences

wo_status_model <- lmer(dN ~ poly(julian_day,2) + (1| whale_id), data=Bp_isotopic_data_use)
status_model <- lmer(dN ~ poly(julian_day,2) + as.factor(Bp_isotopic_data_use$Status) + (1| whale_id), data=Bp_isotopic_data_use)

anova(wo_status_model, status_model, test="Chisq")

### Decision: Evidence that status does not improve the model. Do not include.

#--------------------------------------------------------------------------------------
### Sliding windows
#--------------------------------------------------------------------------------------

### 1) Weather variables 
# 1.1) Relative
#############################################################################################################

library(climwin)

all_env_wea_data1 <- scale(all_env_wea_data[,c("thetao_ice", "thetao_atl", "chl_ice","chl_atl")], center=TRUE, scale=TRUE) #comprovar que s'agafen les variables
all_env_wea_data1 <- as.data.frame(all_env_wea_data1)
all_env_wea_data$thetao_ice <- all_env_wea_data1$thetao_ice
all_env_wea_data$thetao_atl <- all_env_wea_data1$thetao_atl
all_env_wea_data$chl_ice <- all_env_wea_data1$chl_ice
all_env_wea_data$chl_atl <- all_env_wea_data1$chl_atl

Bp_dN_sw_weather_m_rel <- climwin::slidingwin(baseline = lme4::lmer(dN ~ 1 + poly(Bp_isotopic_data_use$julian_day,2) + (1| Bp_isotopic_data_use$whale_id), data=Bp_isotopic_data_use),  ### this is the base model to which the program adds the climate variables
  xvar = list(NAO=all_env_wea_data$NAO_index,
              AMO= all_env_wea_data$AMO_index,
              AMOC= all_env_wea_data$AMOC_index,
              SST_ice=all_env_wea_data$thetao_ice,
              SST_atl= all_env_wea_data$thetao_atl,
              CHL_ice= all_env_wea_data$chl_ice,
              CHL_atl= all_env_wea_data$chl_atl
              ),
  type = "relative", 
  range= c(36,0), 
  stat = c("mean"),
  func = c("lin"),
  cinterval = "month",
  cmissing = "method2", #només afecta a AMOC index on fa la mitjana en els darrers dos mesos
  cdate = all_env_wea_data$date_use, bdate = Bp_isotopic_data_use$measured_date_formatted)


#length(unique(all_env_wea_data$date_use))
#length(all_env_wea_data$date_use)

# Weather relative results
#nitrogen
Bp_dN_sw_weather_m_rel$combos #check all of the slinding windows in the analyses

#AIC difference < 2

#Amb julian day, 48-0 months

a<-head(Bp_dN_sw_weather_m_rel[[1]]$Dataset, 20) ## Lin_func Mean NAO, close window-open window (34-31) (beta coefficient= -0.36)
a <- a %>% filter(deltaAICc < (a$deltaAICc[1])+2)
nrow_a <- nrow(a)
a$Index <- c(rep("NAO", nrow_a))
b<-head(Bp_dN_sw_weather_m_rel[[2]]$Dataset, 20) ## Lin_func Mean AMO, close window-open window (32-27) (beta coefficient= 2.57)
b <- b %>% filter(deltaAICc < (b$deltaAICc[1])+2)
nrow_b <- nrow(b)
b$Index <- c(rep("AMO", nrow_b))
c<-head(Bp_dN_sw_weather_m_rel[[3]]$Dataset, 20) ## Lin_func Mean AMOC, close window-open window (30-13) (beta coefficient= 0.35)
c <- c %>% filter(deltaAICc < (c$deltaAICc[1])+2)
nrow_c <- nrow(c)
c$Index <- c(rep("AMOC", nrow_c))
d<-head(Bp_dN_sw_weather_m_rel[[4]]$Dataset, 20) ## Lin_func Mean SST ice, close window-open window (13-2) (beta coefficient= 0.95)
d <- d %>% filter(deltaAICc < (d$deltaAICc[1])+2)
nrow_d <- nrow(d)
d$Index <- c(rep("SST ice", nrow_d))
e<-head(Bp_dN_sw_weather_m_rel[[5]]$Dataset, 20) ## Lin_func Mean SST atl, close window-open window (48-37) (beta coefficient= -1.76)
e <- e %>% filter(deltaAICc < (e$deltaAICc[1])+2)
nrow_e <- nrow(e)
e$Index <- c(rep("SST atl", nrow_e))
f<-head(Bp_dN_sw_weather_m_rel[[6]]$Dataset, 20) ## Lin_func Mean CHL ice, close window-open window (33-0) (beta coefficient= -11.14)
f <- f %>% filter(deltaAICc < (f$deltaAICc[1])+2)
nrow_f <- nrow(f)
f$Index <- c(rep("CHL ice", nrow_f))
g<-head(Bp_dN_sw_weather_m_rel[[7]]$Dataset, 20) ## Lin_func Mean CHL atl, close window-open window (48-32) (beta coefficient= 16.25)
g <- g %>% filter(deltaAICc < (g$deltaAICc[1])+2)
nrow_g <- nrow(g)
g$Index <- c(rep("CHL atl", nrow_g))

data <- rbind(a,b,c,d,e,f,g)

levels=c("NAO","AMO","AMOC", "SST ice","SST atl","CHL ice","CHL atl")

data$Index <- factor(data$Index, levels=rev(c("NAO","AMO","AMOC", "SST ice","SST atl","CHL ice","CHL atl")))

ggplot(data) +
  geom_vline(aes(xintercept = data$Index), color = "gray", linetype = "dotted") +
  geom_linerange(aes(x = data$Index, ymin = WindowOpen, ymax = WindowClose), color = '#B3CDE3', size = 20, alpha = 0.5) +  # Afegeix les línies de separació
  geom_hline(yintercept = seq(1,48,1), lty = 1, lwd = 0.1, colour = "white") +
  geom_hline(yintercept = 0, lty = 1, lwd = 1, colour = "yellow") +
  coord_flip() +
  labs(x = "", y = "Time from last keratin sample (months)", title = "") +
  theme_classic() + ylim(48,0)+
  scale_x_discrete(breaks=levels,labels=c("Avg NAO index","Avg AMO index","Avg AMOC index","Avg SST (°C)", "Avg SST Atl. (°C)", "Avg CHL (mg m^-3)", "Avg CHL Atl. (mg m^-3)")) + theme(aspect.ratio=3/4)

#--------------------------------------------------------------------------------------
### Randomizations
#--------------------------------------------------------------------------------------

### 1) Weather variables 

###1.1. Relative weather variables
Bp_dN_rand_weather_m <- randwin(repeats = 1000, #100
                                baseline = lme4::lmer(dN ~ 1 + poly(Bp_isotopic_data_use$julian_day,2) + (1| Bp_isotopic_data_use$whale_id), data=Bp_isotopic_data_use),  ### this is the base model to which the program adds the climate variables
                                xvar = list(NAO=all_env_wea_data$NAO_index,
                                            AMO= all_env_wea_data$AMO_index,
                                            AMOC= all_env_wea_data$AMOC_index,
                                            SST_ice=all_env_wea_data$thetao_ice,
                                            SST_atl= all_env_wea_data$thetao_atl,
                                            CHL_ice= all_env_wea_data$chl_ice,
                                            CHL_atl= all_env_wea_data$chl_atl
                                ),
                                type = "relative", 
                                range= c(48,0), 
                                stat = c("mean"),
                                func = c("lin"),
                                cinterval = "month",
                                cmissing = "method2", #només afecta a AMOC index on fa la mitjana en els darrers dos mesos
                                cdate = all_env_wea_data$date_use, bdate = Bp_isotopic_data_use$measured_date_formatted)

Bp_dN_rand_weather_m


#########
### Compare randomizations to original sliding window
#########

# Mean NAO
pvalue(dataset = Bp_dN_sw_weather_m_rel[[1]]$Dataset, datasetrand = Bp_dN_rand_weather_m[[1]], 
       metric="AIC", sample.size=12) # "<0.001"

# Mean AMO
pvalue(dataset = Bp_dN_sw_weather_m_rel[[2]]$Dataset, datasetrand = Bp_dN_rand_weather_m[[2]], 
       metric="AIC", sample.size=12) # "<0.001"

# Mean AMOC
pvalue(dataset = Bp_dN_sw_weather_m_rel[[3]]$Dataset, datasetrand = Bp_dN_rand_weather_m[[3]], 
       metric="AIC", sample.size=12) # "<0.001"

# Mean SST ice
pvalue(dataset = Bp_dN_sw_weather_m_rel[[4]]$Dataset, datasetrand = Bp_dN_rand_weather_m[[4]], 
       metric="AIC", sample.size=12) # "<0.001"

# Mean SST atl
pvalue(dataset = Bp_dN_sw_weather_m_rel[[5]]$Dataset, datasetrand = Bp_dN_rand_weather_m[[5]], 
       metric="AIC", sample.size=12) # "<0.001"

# Mean CHL ice
pvalue(dataset = Bp_dN_sw_weather_m_rel[[6]]$Dataset, datasetrand = Bp_dN_rand_weather_m[[6]], 
       metric="AIC", sample.size=12) # "<0.001"

# Mean CHL atl
pvalue(dataset = Bp_dN_sw_weather_m_rel[[7]]$Dataset, datasetrand = Bp_dN_rand_weather_m[[7]], 
       metric="AIC", sample.size=12) # "<0.001"



###########################################################################################
### Extract weather variables from time periods identified by sliding window ########
###########################################################################################
### Weather

###############
### NAO
###############

### Bp dN: Best window for mean NAO linear 34-31

Bp_isotopic_data_use$mean_NAO_linear_start_m <- ymd(Bp_isotopic_data_use$measured)%m-% months(34)
Bp_isotopic_data_use$mean_NAO_linear_end_m <- ymd(Bp_isotopic_data_use$measured)%m-% months(31)


###############
### AMO
###############

### Bp dN: Best window for mean AO linear 52-35

Bp_isotopic_data_use$mean_AMO_linear_start_m <- ymd(Bp_isotopic_data_use$measured)%m-% months(32)
Bp_isotopic_data_use$mean_AMO_linear_end_m <- ymd(Bp_isotopic_data_use$measured)%m-% months(27)


###############
### AMOC
###############

### Bp dN: Best window for mean AMOC linear 42-30

Bp_isotopic_data_use$mean_AMOC_linear_start_m <- ymd(Bp_isotopic_data_use$measured)%m-% months(30)
Bp_isotopic_data_use$mean_AMOC_linear_end_m <- ymd(Bp_isotopic_data_use$measured)%m-% months(13)


###############
### SST ice
###############

### Bp dN: Best window for mean CHL linear 33-3

Bp_isotopic_data_use$mean_SST_ice_linear_start_m <- ymd(Bp_isotopic_data_use$measured)%m-% months(11)
Bp_isotopic_data_use$mean_SST_ice_linear_end_m <- ymd(Bp_isotopic_data_use$measured)%m-% months(1)


###############
### SST atl
###############

### Bp dN: Best window for mean THETAO linear 53-32

Bp_isotopic_data_use$mean_SST_atl_linear_start_m <- ymd(Bp_isotopic_data_use$measured)%m-% months(34)
Bp_isotopic_data_use$mean_SST_atl_linear_end_m <- ymd(Bp_isotopic_data_use$measured)%m-% months(0)


###############
### CHL ice
###############

### Bp dN: Best window for mean CHL_Atl linear 48-34

Bp_isotopic_data_use$mean_CHL_ice_linear_start_m <- ymd(Bp_isotopic_data_use$measured)%m-% months(32)
Bp_isotopic_data_use$mean_CHL_ice_linear_end_m <- ymd(Bp_isotopic_data_use$measured)%m-% months(4)


###############
### CHL atl
###############

### Bp dN: Best window for mean THETAO linear 48-37

Bp_isotopic_data_use$mean_CHL_atl_linear_start_m <- ymd(Bp_isotopic_data_use$measured)%m-% months(27)
Bp_isotopic_data_use$mean_CHL_atl_linear_end_m <- ymd(Bp_isotopic_data_use$measured)%m-% months(16)


###############
### Select entire 30 day period for all weather variables for summary table
###############

Bp_isotopic_data_use$overall_start <- ymd(Bp_isotopic_data_use$measured)%m-% months(36)
Bp_isotopic_data_use$overall_end <- ymd(Bp_isotopic_data_use$measured)%m-% months(0)

###########################################
### Convert dates to proper posix format
###########################################

### Convert weather dates to POSIXct or else the following code will not work 
### (i.e., it won't work with same date format as sliding window)

all_env_wea_data$date.int <- strptime(paste(all_env_wea_data$date), format="%Y-%m-%d") 
all_env_wea_data$date_new <- as.POSIXct(all_env_wea_data$date.int)

head(all_env_wea_data)

### Remove date.int or else loop below will not work

all_env_wea_data <- all_env_wea_data[,-18] #assegurar-se que s'esborra la columna date_int
all_env_wea_data

#---------------------------------------------------------------------------------
#Nao
#---------------------------------------------------------------------------------
###############
### Mean NAO
###############

Bp_isotopic_data_use$mean_NAO_start_m_int <- strptime(paste(Bp_isotopic_data_use$mean_NAO_linear_start_m), format="%Y-%m-%d") 
Bp_isotopic_data_use$mean_NAO_start_m_date <- as.POSIXct(Bp_isotopic_data_use$mean_NAO_start_m_int)

Bp_isotopic_data_use$mean_NAO_end_m_int <- strptime(paste(Bp_isotopic_data_use$mean_NAO_linear_end_m), format="%Y-%m-%d") 
Bp_isotopic_data_use$mean_NAO_end_m_date <- as.POSIXct(Bp_isotopic_data_use$mean_NAO_end_m_int)


#---------------------------------------------------------------------------------
#Amo
#---------------------------------------------------------------------------------
###############
### Mean AMO
###############

Bp_isotopic_data_use$mean_AMO_start_m_int <- strptime(paste(Bp_isotopic_data_use$mean_AMO_linear_start_m), format="%Y-%m-%d") 
Bp_isotopic_data_use$mean_AMO_start_m_date <- as.POSIXct(Bp_isotopic_data_use$mean_AMO_start_m_int)

Bp_isotopic_data_use$mean_AMO_end_m_int <- strptime(paste(Bp_isotopic_data_use$mean_AMO_linear_end_m), format="%Y-%m-%d") 
Bp_isotopic_data_use$mean_AMO_end_m_date <- as.POSIXct(Bp_isotopic_data_use$mean_AMO_end_m_int)


#---------------------------------------------------------------------------------
#Amoc
#---------------------------------------------------------------------------------
###############
### Mean AMOC
###############

Bp_isotopic_data_use$mean_AMOC_start_m_int <- strptime(paste(Bp_isotopic_data_use$mean_AMOC_linear_start_m), format="%Y-%m-%d") 
Bp_isotopic_data_use$mean_AMOC_start_m_date <- as.POSIXct(Bp_isotopic_data_use$mean_AMOC_start_m_int)

Bp_isotopic_data_use$mean_AMOC_end_m_int <- strptime(paste(Bp_isotopic_data_use$mean_AMOC_linear_end_m), format="%Y-%m-%d") 
Bp_isotopic_data_use$mean_AMOC_end_m_date <- as.POSIXct(Bp_isotopic_data_use$mean_AMOC_end_m_int)


#---------------------------------------------------------------------------------
#SST ice
#---------------------------------------------------------------------------------
###############
### Mean SST ice
###############

Bp_isotopic_data_use$mean_SST_ice_start_m_int <- strptime(paste(Bp_isotopic_data_use$mean_SST_ice_linear_start_m), format="%Y-%m-%d") 
Bp_isotopic_data_use$mean_SST_ice_start_m_date <- as.POSIXct(Bp_isotopic_data_use$mean_SST_ice_start_m_int)

Bp_isotopic_data_use$mean_SST_ice_end_m_int <- strptime(paste(Bp_isotopic_data_use$mean_SST_ice_linear_end_m), format="%Y-%m-%d") 
Bp_isotopic_data_use$mean_SST_ice_end_m_date <- as.POSIXct(Bp_isotopic_data_use$mean_SST_ice_end_m_int)


#---------------------------------------------------------------------------------
#SST atl
#---------------------------------------------------------------------------------
###############
### Mean SST atl
###############

Bp_isotopic_data_use$mean_SST_atl_start_m_int <- strptime(paste(Bp_isotopic_data_use$mean_SST_atl_linear_start_m), format="%Y-%m-%d") 
Bp_isotopic_data_use$mean_SST_atl_start_m_date <- as.POSIXct(Bp_isotopic_data_use$mean_SST_atl_start_m_int)

Bp_isotopic_data_use$mean_SST_atl_end_m_int <- strptime(paste(Bp_isotopic_data_use$mean_SST_atl_linear_end_m), format="%Y-%m-%d") 
Bp_isotopic_data_use$mean_SST_atl_end_m_date <- as.POSIXct(Bp_isotopic_data_use$mean_SST_atl_end_m_int)


#---------------------------------------------------------------------------------
#CHL ice
#---------------------------------------------------------------------------------
###############
### Mean CHL ice
###############

Bp_isotopic_data_use$mean_CHL_ice_start_m_int <- strptime(paste(Bp_isotopic_data_use$mean_CHL_ice_linear_start_m), format="%Y-%m-%d") 
Bp_isotopic_data_use$mean_CHL_ice_start_m_date <- as.POSIXct(Bp_isotopic_data_use$mean_CHL_ice_start_m_int)

Bp_isotopic_data_use$mean_CHL_ice_end_m_int <- strptime(paste(Bp_isotopic_data_use$mean_CHL_ice_linear_end_m), format="%Y-%m-%d") 
Bp_isotopic_data_use$mean_CHL_ice_end_m_date <- as.POSIXct(Bp_isotopic_data_use$mean_CHL_ice_end_m_int)


#---------------------------------------------------------------------------------
#CHL atl
#---------------------------------------------------------------------------------
###############
### Mean CHL atl
###############

Bp_isotopic_data_use$mean_CHL_atl_start_m_int <- strptime(paste(Bp_isotopic_data_use$mean_CHL_atl_linear_start_m), format="%Y-%m-%d") 
Bp_isotopic_data_use$mean_CHL_atl_start_m_date <- as.POSIXct(Bp_isotopic_data_use$mean_CHL_atl_start_m_int)

Bp_isotopic_data_use$mean_CHL_atl_end_m_int <- strptime(paste(Bp_isotopic_data_use$mean_CHL_atl_linear_end_m), format="%Y-%m-%d") 
Bp_isotopic_data_use$mean_CHL_atl_end_m_date <- as.POSIXct(Bp_isotopic_data_use$mean_CHL_atl_end_m_int)


#---------------------------------------------------------------------------------

### Overal time period

Bp_isotopic_data_use$overall_start_int <- strptime(paste(Bp_isotopic_data_use$overall_start), format="%Y-%m-%d") 
Bp_isotopic_data_use$overall_start_date <- as.POSIXct(Bp_isotopic_data_use$overall_start_int)

Bp_isotopic_data_use$overall_end_int <- strptime(paste(Bp_isotopic_data_use$overall_end), format="%Y-%m-%d") 
Bp_isotopic_data_use$overall_end_date <- as.POSIXct(Bp_isotopic_data_use$overall_end_int)



#############################################
### Create loop to extract weather variables for each sample and save as a csv
#############################################

### Empty lists

mean_NAO_dN_list <- list()


mean_AMO_dN_list <- list()


mean_AMOC_dN_list <- list()


mean_SST_ice_dN_list <- list()


mean_SST_atl_dN_list <- list()


mean_CHL_ice_dN_list <- list()


mean_CHL_atl_dN_list <- list()


overall_list <- list()


for (i in 1:1383) { #numero de samples 
  #mean NAO
  mean_NAO_dN_list[[i]] <- all_env_wea_data %>% 
    filter(date_new >= Bp_isotopic_data_use$mean_NAO_start_m_date[i] 
           & date_new <= Bp_isotopic_data_use$mean_NAO_end_m_date[i]) %>% 
    
    summarise(mean_NAO_dN = mean(NAO_index))
  
  
  #mean AMO
  mean_AMO_dN_list[[i]] <- all_env_wea_data %>% 
    filter(date_new >= Bp_isotopic_data_use$mean_AMO_start_m_date[i] 
           & date_new <= Bp_isotopic_data_use$mean_AMO_end_m_date[i]) %>% 
    
    summarise(mean_AMO_dN = mean(AMO_index))
  
  
  
  #mean AMOC
  mean_AMOC_dN_list[[i]] <- all_env_wea_data %>% 
    filter(date_new >= Bp_isotopic_data_use$mean_AMOC_start_m_date[i] 
           & date_new <= Bp_isotopic_data_use$mean_AMOC_end_m_date[i]) %>% 
    
    summarise(mean_AMOC_dN = mean(AMOC_index))
  
  
  #mean SST ice
  mean_SST_ice_dN_list[[i]] <- all_env_wea_data %>% 
    dplyr::filter(date_new >= Bp_isotopic_data_use$mean_SST_ice_start_m_date[i] 
                  & date_new <= Bp_isotopic_data_use$mean_SST_ice_end_m_date[i]) %>% 
    
    dplyr::summarise(mean_SST_ice_dN = mean(thetao_ice, na.rm = TRUE))
  
  
  #mean SST atl
  mean_SST_atl_dN_list[[i]] <- all_env_wea_data %>% 
    filter(date_new >= Bp_isotopic_data_use$mean_SST_atl_start_m_date[i] 
           & date_new <= Bp_isotopic_data_use$mean_SST_atl_end_m_date[i]) %>% 
    
    summarise(mean_SST_atl_dN = mean(thetao_atl))
  
  
  #mean CHL ice
  mean_CHL_ice_dN_list[[i]] <- all_env_wea_data %>% 
    dplyr::filter(date_new >= Bp_isotopic_data_use$mean_CHL_ice_start_m_date[i] 
                  & date_new <= Bp_isotopic_data_use$mean_CHL_ice_end_m_date[i]) %>% 
    
    dplyr::summarise(mean_CHL_ice_dN = mean(chl_ice, na.rm = TRUE))
  
  
  
  #mean CHL atl
  mean_CHL_atl_dN_list[[i]] <- all_env_wea_data %>% 
    filter(date_new >= Bp_isotopic_data_use$mean_CHL_atl_start_m_date[i] 
           & date_new <= Bp_isotopic_data_use$mean_CHL_atl_end_m_date[i]) %>% 
    
    summarise(mean_CHL_atl_dN = mean(chl_atl))
  
  overall_list[[i]] <- all_env_wea_data %>% 
    filter(date_new >= Bp_isotopic_data_use$overall_start_date[i] 
           & date_new <= Bp_isotopic_data_use$overall_end_date[i]) %>% 
    
    summarise(mean_NAO = mean(NAO_index, na.rm = TRUE),
              mean_AMO = mean(AMO_index, na.rm = TRUE),
              mean_AMOC = mean(ENSO3_4_index, na.rm = TRUE),
              mean_SST_ice = mean(thetao_ice, na.rm = TRUE),
              mean_SST_atl = mean(thetao_atl, na.rm = TRUE),
              mean_CHL_ice = mean(chl_ice, na.rm = TRUE),
              mean_CHL_atl = mean(chl_atl, na.rm = TRUE))
  
  
  ### Convert lists to rows per sample ID
  
  mean_NAO_dN_data <- bind_rows(mean_NAO_dN_list)
  
  mean_AMO_dN_data <- bind_rows(mean_AMO_dN_list)
  
  mean_AMOC_dN_data <- bind_rows(mean_AMOC_dN_list)
  
  mean_SST_ice_dN_data <- bind_rows(mean_SST_ice_dN_list)
  
  mean_SST_atl_dN_data <- bind_rows(mean_SST_atl_dN_list)
  
  mean_CHL_ice_dN_data <- bind_rows(mean_CHL_ice_dN_list)
  
  mean_CHL_atl_dN_data <- bind_rows(mean_CHL_atl_dN_list)
  
  overall_data <- bind_rows(overall_list)
  
  
  ### Convert each variable to a column in a larger dataset
  Bp_dN_windows <- bind_cols(mean_NAO_dN_data,
                             mean_AMO_dN_data, 
                             mean_AMOC_dN_data, 
                             mean_SST_ice_dN_data, 
                             mean_SST_atl_dN_data,
                             mean_CHL_ice_dN_data,
                             mean_CHL_atl_dN_data,
                             overall_data)
  
}

head(Bp_dN_windows)


### Combine with nestquadg measurement data

head(Bp_isotopic_data_use)

Bp_isotopic_data_new <- Bp_isotopic_data_use # Extract only columns of interest from original data

Bp_dN_final <- bind_cols(Bp_isotopic_data_new, Bp_dN_windows)


### Save data

write.csv(Bp_dN_final, "/Users/marcruizisagales/Desktop/Doctorat/Analisis_isotops_barbes/Projecte_barbes_clima/Environmental_vars/Bp_isotopic_weather_summarized_final_nitrogen_13_nov_mean_seasonality.csv")

