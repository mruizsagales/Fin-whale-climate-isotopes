#--------------------------------------------------------------------------------
#Figure 2 - Example of time-estimation in baleen plates
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
  legend.title = element_text(size = 14), # Legend title
  legend.text = element_text(size = 12),  # Legend text
  #strip.text = element_text(size = 15),    # Facet labels
  aspect.ratio = 3/4
)

#2. Import data

df <- read_excel("/Users/marcruizisagales/Documents/GitHub/Climate-baleen-plate-isotopes/data/0_point_alignments_29_juny_total_-4.xlsx") # Import dataset of the one-centimetre-spaced stable isotope data along the baleen plate of fin whales

min_cm_minus_4 <- whale_isos %>% dplyr::group_by(Whale) %>% dplyr::summarise(min_cm = min(Cm)) # Make sure that all the stable isotope data start at the position -4 cm 

fin_whale <- filter(df, Whale == "F22046")

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# convert sample number to days based on growth rate of baleen
growth_rate <- 16  # baleen growth rate in centimetres per year
last_sample <- ymd("2022-07-22") # final date of last keratin sample

fin_whale <- 
  fin_whale %>%
  mutate(days = Cm * 365.25 / growth_rate) %>%
  mutate(days = days - days[1]) %>%
  mutate(rev_days = days - days[length(days)]) %>%
  mutate(sample_date = last_sample + rev_days) %>%
  mutate(year = year(sample_date))

# Create a variable to label the years
year_labels <- with(fin_whale, seq(min(year), max(year), by = 1))
year_breaks <- which(duplicated(fin_whale$year) == FALSE)
year_labels <- c(2022, 2021,2020)
year_breaks <- c(45-9,45-25,45-41)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# plot the isotope data
point_size <- 5

isop <- ggplot(fin_whale, aes(x = rev(Cm))) + 
  theme_classic(base_size = 16) + 
  xlab("Baleen segment from gingiva (Cm)") + 
  scale_x_continuous(breaks = seq(0-4, 50-4, 10), 
                     labels = seq(50-4, 0-4, -10),
                     sec.axis = sec_axis(~., 
                                         breaks = year_breaks,
                                         labels = year_labels,
                                         name = "Year"))

isop <- isop + geom_line(aes(y = dC), color="black") + 
  geom_point(aes(y = dC), size = point_size, shape = 21, fill = "black") +
  ylab(expression(paste(delta^{13}, "C (\u2030)")))

# add the d15N data but need to scale it to d13C
a <- mean(na.omit(fin_whale$dN))# add the d15N data but need to scale it to d13C
a
b <- mean(na.omit(fin_whale$dC))# add the d15N data but need to scale it to d13C
b
v_shift <- a - b

isop <- isop + geom_line(aes(y = dN - v_shift), color="grey") + 
  geom_point(aes(y = (dN - v_shift)), size = point_size, 
             shape = 21, fill = "grey")

# add the second axis and rescale it appropriately
isop <- isop + 
  scale_y_continuous(
    sec.axis = sec_axis(~.+v_shift, 
                        name = expression(paste(delta^{15}, "N (\u2030)"))))

isop <- isop + geom_vline(xintercept = rev(year_breaks),
                          color = "grey", lty = 2) + theme_article(base_size = 15) + standard_theme 

print(isop) 

ggsave("/Users/marcruizisagales/Documents/GitHub/Climate-baleen-plate-isotopes/png/Figure_2_Example_time-estimation_baleen_plates.png", last_plot(), 
       dpi = 300,  width = 16, height = 12, units = "cm")
       

