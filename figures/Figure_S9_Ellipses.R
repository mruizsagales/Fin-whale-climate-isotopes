
#--------------------------------------------------------------------------------
# Figure S9 - Ellipses and SEAc trough time
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
library(devtools)
library(SIBER)
library(egg)

# 2. Import stable isotope data
df <- read_excel("~/Desktop/Doctorat/Analisis_isotops_barbes/Projecte_barbes_clima/Environmental_vars/All_merged_time_calibrated_2013_to_2022_Suess_cor.xlsx") # Import Suess-corrected dataset of stable isotope data along the baleen plate of fin whales 

df_1 <- df[complete.cases(df$dN),] # remove NA
df_2 <- df_1[complete.cases(df_1$d13c.cor),] # remove NA

sort(unique(df_2$Year_from_sample_date)) # years


#df_2 <- subset(df_2, !(Year_from_sample_date == 2010 & Whale == "F13065") & !(Year_from_sample_date == 2010 & Whale == "F13076")& !(Year_from_sample_date == 2012 & Whale == "F15083")& !(Year_from_sample_date == 2012 & Whale == "F15097")) # remove

# 1. Ellipses

df_3 <- df_2[, c("d13c.cor", "dN", "Year_from_sample_date", "Species")]
names(df_3) <- c("iso1","iso2", "group", "community")
df_4 <- as.data.frame(df_3)
# create the siber object
#siber.example <- createSiberObject(demo.siber.data)

siber.fin <- createSiberObject(df_4)

# Calculate summary statistics for each group: TA, SEA and SEAc
group.ML <- groupMetricsML(siber.fin)
print(group.ML)
df_group.ML <- as.data.frame(group.ML)
df_group.ML$RowNames <- rownames(df_group.ML)
# Reshape the dataframe
reshaped_dataframe <- gather(df_group.ML, key = "Key", value = "Value", -RowNames)

# Separate the 'Key' column into 'Whale' and 'Year'
reshaped_dataframe <- separate(reshaped_dataframe, Key, into = c("Whale", "Year"), sep = "\\.", remove = FALSE)

# Convert the dataframe to the desired format
library(tidyverse)

new_df <- as.data.frame(group.ML)
names(new_df) <- c(2010:2022)
# Transpose the dataframe and convert to tidy format
tidy_df <- t(new_df) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Year") %>%
  gather(key = "Variable", value = "Value", -Year)

# Print the tidy dataframe
print(tidy_df)

# Communuty metrics ML
community.ML <- communityMetricsML(siber.fin) 
print(community.ML)

###
sbg_1 <- df %>% 
  dplyr::group_by(Year_from_sample_date, Whale) %>% 
  dplyr::summarise(mC = mean(d13c.cor, na.rm=TRUE), 
                   sdC = sd(d13c.cor, na.rm=TRUE), 
                   mN = mean(dN, na.rm=TRUE), 
                   sdN = sd(dN, na.rm=TRUE) )


# set colors
coul <- brewer.pal(9, "RdYlBu") 
coul <- colorRampPalette(coul)(30)

# plot points

values <- setNames(coul, sort(unique(sbg_1$Whale)))
first.plot <- ggplot(data = as.data.frame(df), 
                     mapping = aes(x = d13c.cor, 
                                   y = dN)) + 
  geom_point(aes(color = Whale), size = 0.22, alpha=0.9) +
  facet_wrap(~Year_from_sample_date, ncol = 4) +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) + 
  theme(text = element_text(size=20)) + theme_article(base_size = 16) +
  scale_colour_manual(values = values)

print(first.plot)


# add errorbars
second.plot <- first.plot + 
  geom_errorbar(data = sbg_1, mapping = aes(x = mC, y = mN,
                                             ymin = mN - 1.96*sdN, 
                                             ymax = mN + 1.96*sdN, col = sbg_1$Whale),width=0.2) +
  geom_errorbar(data = sbg_1, mapping = aes(x = mC, y = mN,
                                            xmin = mC - 1.96*sdC,
                                            xmax = mC + 1.96*sdC, col = sbg_1$Whale), width = 0.2) + 
  geom_point(data = sbg_1, aes(x = mC, 
                             y = mN, col = sbg_1$Whale), shape = 22, size = 1,
             alpha = 0.7, show.legend = FALSE) + coord_fixed(ratio=3/8) +
  scale_colour_manual(values = values)
second.plot

# addellipses

p.ell <- 0.95 # decide how big an ellipse you want to draw

ellipse.plot <- second.plot + 
  stat_ellipse(aes(group = Year_from_sample_date), alpha=0.01, 
               linetype = 2,
               color="darkgrey",
               level = p.ell,
               type = "norm",
               geom = "polygon")+ ylim(6,14)

print(ellipse.plot) 

print(ellipse.plot) + theme(legend.position = c(0.6, 0.1),
                            legend.direction = "horizontal", legend.title = element_blank()) + ### Set aspect ratio of the graph n/n = square
  theme(aspect.ratio=2/4) 


# 2. Plot 2

sbg_2 <- df_2 %>% 
  dplyr::group_by(Year_from_sample_date) %>% 
  dplyr::summarise(count = n(),
                   mC = mean(d13c.cor), 
                   sdC = sd(d13c.cor), 
                   mN = mean(dN), 
                   sdN = sd(dN))

colnames(sbg_2) <- c("Year","count","mC","sdC","mN","sdN")
tidy_df_SEAc <- filter(tidy_df, Variable == "SEAc")

Seac <- merge(sbg_2, tidy_df_SEAc, by= "Year")



Seac$Year <- as.double(Seac$Year)
Seac$Value <- as.numeric(Seac$Value)
Seac$count <- as.numeric(Seac$count)

Seac_gam <- mgcv::gam(Value ~  s(Year),
                        data= Seac)

plot(Value ~ Year, Seac)

summary(Seac_gam)

par(mar=c(5,5,5,5))
plot.gam(Seac_gam)

model_2_p <- predict_gam(SDn_mN_gam)
model_2_p

model_2_p %>%
  ggplot(aes(mN, count, z = fit)) +
  geom_raster(aes(fill = fit)) +
  geom_contour(colour = "white") +
  scale_fill_continuous(name = "y") +
  theme_minimal() +
  theme(legend.position = "top")

predict_gam(SDn_mN_gam, exclude_terms = "s(count)") %>%
  ggplot(aes(mN, fit)) + geom_smooth_ci() + geom_line(color = "darkgrey", linewidth=2) +
  geom_ribbon(aes(ymin = fit-1.96 *se.fit, ymax = fit+1.96 *se.fit), alpha = 0,color = "black", linetype = "dotted") + 
  ylab(expression(paste("Annual ",delta^{15}, "N (\u2030) SD"))) +
  xlab(expression(paste("Annual ",delta^{15}, "N (\u2030) mean"))) + theme_article(base_size = 20)

# 3. Plot 3

sbg_2$Year_from_sample_date <- as.numeric(sbg_2$Year_from_sample_date)
sbg_2$sdN <- as.numeric(sbg_2$sdN)
sbg_2$mN <- as.numeric(sbg_2$mN)

SDn_mN_gam <- mgcv::gam(sdN ~  s(mN) + s(count),
                                data= sbg_2,method="REML", select=TRUE)

summary(SDn_mN_gam)
par(mar=c(5,5,5,5))
plot.gam(SDn_mN_gam)

model_2_p <- predict_gam(SDn_mN_gam)
model_2_p

model_2_p %>%
  ggplot(aes(mN, count, z = fit)) +
  geom_raster(aes(fill = fit)) +
  geom_contour(colour = "white") +
  scale_fill_continuous(name = "y") +
  theme_minimal() +
  theme(legend.position = "top")

predict_gam(SDn_mN_gam, exclude_terms = "s(count)") %>%
  ggplot(aes(mN, fit)) + geom_smooth_ci() + geom_line(color = "darkgrey", linewidth=2) +
  geom_ribbon(aes(ymin = fit-1.96 *se.fit, ymax = fit+1.96 *se.fit), alpha = 0,color = "black", linetype = "dotted") + 
  ylab(expression(paste("Annual ",delta^{15}, "N (\u2030) SD"))) +
  xlab(expression(paste("Annual ",delta^{15}, "N (\u2030) mean"))) + theme_article(base_size = 20)

#####

