#--------------------------------------------------------------------------------
#Oscilacions de la senyal isotòpica en la barba de rorqual comú (Islàndia 2010-2022)
#--------------------------------------------------------------------------------

# Inspirat en script de Clive Trueman
# Adaptat per: Andrew Jackson June 2017 & Natalie Cooper Oct 2017.
# Readaptat per: Marc Ruiz-Sagalés

#--------------------------------------------------------------------------------

#IMPORTANT: S'ha de treure els hibrids i immadurs? per al paper de les barbes

# Load libraries
library(tidyverse)
library(lubridate)
library(ggrepel)
library(readxl)
library(ggplot2)
library(egg)
library(zoo)
library(TSA)

#Import data

#whale_isos <-read_excel("~/Desktop/Doctorat/Analisis_isotops_barbes/Clima/Data/Dades_barbes.xlsx", 
                        col_types = c("numeric", "numeric", "numeric", 
                                      "text", "text", "text", "text", "numeric", 
                                      "text", "numeric", "numeric", "date", 
                                      "text", "text", "text", "numeric", 
                                      "numeric", "numeric", "text", "numeric", 
                                      "numeric", "numeric", "numeric", 
                                      "numeric", "numeric", "numeric"))

#dataset amb el 0 ben alienat en totes les barbes
whale_isos <- read_excel("~/Desktop/Doctorat/Analisis_isotops_barbes/0_point_alignments_29_juny_total_-4.xlsx")

# Author: Clive Trueman
# Adapted by: Andrew Jackson June 2017 & Natalie Cooper Oct 2017.
# About: Script to generate figure 1

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Load libraries
library(tidyverse)
library(lubridate)
library(ggrepel)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Figure 1 - raw whale isotope data with two y-axes
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# read in data and extract just the blue whale

KC7 <- filter(whale_isos, Whale == "F22046")

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# convert sample number to days based on growth rate of baleen
growth_rate <- 16  # baleen growth rate in centimetres per year
last_sample <- ymd("2022-07-22") # final date of last keratin sample

KC7 <- 
  KC7 %>%
  mutate(days = Cm * 365.25 / growth_rate) %>%
  mutate(days = days - days[1]) %>%
  mutate(rev_days = days - days[length(days)]) %>%
  mutate(sample_date = last_sample + rev_days) %>%
  mutate(year = year(sample_date))

# Create a variable to label the years
year_labels <- with(KC7, seq(min(year), max(year), by = 1))
year_breaks <- which(duplicated(KC7$year) == FALSE)
year_labels <- c(2022, 2021,2020)
year_breaks <- c(45-9,45-25,45-41)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# plot the isotope data
point_size <- 4

isop <- ggplot(KC7, aes(x = rev(Cm))) + 
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
a <- mean(na.omit(KC7$dN))# add the d15N data but need to scale it to d13C
a
b <- mean(na.omit(KC7$dC))# add the d15N data but need to scale it to d13C
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
                          color = "grey", lty = 2)

print(isop)

ggsave("manuscript/PeerJ/figures/Figure-1-raw-dC-dN-data.png", isop, 
       device = png(width = 600, height = 400))


# Load necessary library
library(stats)
#NITROGEN
# Get unique whale names
whale_names <- unique(df$Whale)

# Create an empty matrix to store max correlation values
num_whales <- length(whale_names)
correlation_matrix <- matrix(1,nrow = num_whales, ncol = num_whales, dimnames = list(whale_names, whale_names))

# Loop through all possible pairs of whales

for (i in 1:num_whales) {
  for (j in 1:num_whales) {
    if (i != j) {
      # Extract dN for the current pair of whales
      dN_i <- df[df$Whale == whale_names[i], "dN"]
      dN_j <- df[df$Whale == whale_names[j], "dN"]
      # Calculate cross-correlation using ccf function
      ccf_result <- ccf(na.approx(dN_i$dN), na.approx(dN_j$dN), lag.max = 3, plot = FALSE)
      
      # Store the maximum correlation value in the matrix
      correlation_matrix[i, j] <- max(ccf_result$acf)
    }
  }
}

# Print or use the correlation matrix as needed
print(correlation_matrix)
library(corrplot)
corrplot(correlation_matrix, method="color", tl.col="black")
mean(correlation_matrix)

#CARBON
# Get unique whale names
whale_names <- unique(df$Whale)

# Create an empty matrix to store max correlation values
num_whales <- length(whale_names)
correlation_matrix <- matrix(1,nrow = num_whales, ncol = num_whales, dimnames = list(whale_names, whale_names))

# Loop through all possible pairs of whales
for (i in 1:num_whales) {
  for (j in 1:num_whales) {
    if (i != j) {
      # Extract dC for the current pair of whales
      dC_i <- df[df$Whale == whale_names[i], "dC"]
      dC_j <- df[df$Whale == whale_names[j], "dC"]
      # Calculate cross-correlation using ccf function
      ccf_result <- ccf(na.approx(dC_i$dC), na.approx(dC_j$dC), lag.max = 3, plot = FALSE)
      
      # Store the maximum correlation value in the matrix
      correlation_matrix[i, j] <- max(ccf_result$acf)
    }
  }
}

# Print or use the correlation matrix as needed
print(correlation_matrix)
library(corrplot)
corrplot(correlation_matrix, method="color", tl.col="black")
mean(correlation_matrix)

#N vs C
# Get unique whale names
whale_names <- unique(df$Whale)

# Create an empty matrix to store max correlation values
num_whales <- length(whale_names)
correlation_matrix <- matrix(NA,nrow = num_whales, dimnames = list(whale_names))

# Loop through all possible pairs of whales
for (i in 1:num_whales) {
  for (j in 1:num_whales) {
    if (i != j) {
      # Extract dC for the current pair of whales
      dN_i <- df[df$Whale == whale_names[i], "dN"]
      dC_i <- df[df$Whale == whale_names[i], "dC"]
      # Calculate cross-correlation using ccf function
      ccf_result <- ccf(na.approx(dN_i$dN), na.approx(dC_i$dC), lag.max = 3, plot = FALSE)
      
      # Store the maximum correlation value in the matrix
      correlation_matrix[i] <- max(ccf_result$acf)
    }
  }
}

# Print or use the correlation matrix as needed
print(correlation_matrix)
library(corrplot)
color_palette <- colorRampPalette(viridis::viridis(100))(29)
corrplot(correlation_matrix, method="color", tl.col="black", col = rev(color_palette))
mean(correlation_matrix)

whales <- rownames(correlation_matrix)
ccf <- as.vector(correlation_matrix)
ccf_whales <- data.frame(whales, ccf)
range(ccf_whales$ccf)
library(viridis)
ggplot(ccf_whales, aes(x = whales, y = ccf, fill = ccf)) +
  geom_bar(stat = "identity") +
  scale_fill_distiller(palette = "RdBu", direction = 1) +
  scale_y_continuous(limits=c(0,1)) +
  labs(x = "Whales", y="") +
  theme_article(base_size = 15) +
  theme(axis.text.x=element_text(angle=90, hjust=1), legend.position = c(0.9,0.8)) + labs(fill = "CCF")

#________________________________________________________________________________
#Mostres de l'any 2013
#________________________________________________________________________________

whale_isos_2013 <- dplyr::filter(whale_isos, Year == 2013)
unique(whale_isos_2013$Whale) #(n=7)

point_size <- 2
combined_df_2013 <- NULL;

for (i in 1:length(unique(whale_isos_2013$Whale))) {
  df_i <- dplyr::filter(whale_isos, Whale == unique(whale_isos_2013$Whale)[i]) #Filtrem l'individu que volem
  growth_rate <- 16  # baleen growth rate in centimetres per year
  last_sample <- unique(ymd(df_i$Data_capt)) # final date of last keratin sample
  df_i <- 
    df_i %>%
    mutate(days = Cm * 365.25 / growth_rate) %>%
    mutate(days = days - days[1]) %>%
    mutate(rev_days = days - days[length(days)]) %>%
    mutate(sample_date = last_sample + rev_days) %>%
    mutate(year = sample_date)%>%
    mutate(year_rev = rev(sample_date)) #Funció que retrospectivament data cada punt analitzat
  
  df_i$Year_from_sample_date <- format(as.Date(df_i$year_rev), "%Y") #extreure l'any de la data estimada
  df_i$Year_month <- format(as.Date(df_i$year_rev), "%Y/%m") #extreure l'any i el mes de la data estimada
  df_i$Month_from_sample_date <- format(as.Date(df_i$year_rev), "%m") #extreure només el mes de la data estimada
  isop <- ggplot(df_i, aes(x = rev(Cm))) + 
    theme_classic(base_size = 16) + 
    xlab("Baleen sample (cm from youngest sample)") + 
    scale_x_continuous(breaks = df_i$Cm, 
                       labels = rev(df_i$Cm),
                       sec.axis = sec_axis(~., 
                                           breaks = rev(df_i$Cm),
                                           labels = rev(df_i$sample_date),
                                           name = "Date"))
  #Aquest pas es pot saltar si es vol obtenir un plot amb doble y-axis
  isop <- ggplot(df_i, aes(x = rev(sample_date))) + 
    xlab("Date") + 
    scale_x_continuous(breaks = seq(min(df_i$sample_date), max(df_i$sample_date), by = "1 month"),
                       labels = format(seq(min(rev(df_i$sample_date)), max(rev(df_i$sample_date)), by = "1 month"), "%b-%Y"))
  isop <- isop + geom_line(aes(y = dC), col = "black") + 
    geom_point(aes(y = dC), size = point_size, shape = 21, fill = "black", col = "black") +
    ylab(expression(paste(delta^{13}, "C (\u2030)")))+ theme_article()  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  a <- mean(na.omit(df_i$dN))# add the d15N data but need to scale it to d13C
  a
  b <- mean(na.omit(df_i$dC))# add the d15N data but need to scale it to d13C
  b
  v_shift <- a - b
  isop <- isop + geom_line(aes(y = dN - v_shift), col = "darkgrey") + 
    geom_point(aes(y = (dN - v_shift)), size = point_size, 
               shape = 21, fill = "grey", col = "darkgrey")
  isop <- isop + 
    scale_y_continuous(sec.axis = sec_axis(~.+v_shift, 
                                           name = expression(paste(delta^{15}, "N (\u2030)"))),limits = c(-30, -10)) + labs(title= unique(whale_isos_2013$Whale)[i])
  new_df <- df_i
  combined_df_2013 <- rbind(combined_df_2013, new_df)
  print(isop)
  }

library(RColorBrewer)
brewer.pal(n = 6, name = "Reds")
  #Relació de la senyal isotòpica de nitrogen amb la data estimada de cada punt
ggplot() + geom_line(data= combined_df_2013[combined_df_2013$Whale == "F13065",], aes(x= rev(sample_date), y= dN),col="#FEE5D9")+ 
  geom_point(data= combined_df_2013[combined_df_2013$Whale == "F13065",],aes(x= rev(sample_date), y= dN), size = point_size, shape = 21, fill = "#FEE5D9") +
  geom_line(data= combined_df_2013[combined_df_2013$Whale == "F13066",], aes(x= rev(sample_date), y= dN),col="#FCBBA1")+ 
  geom_point(data= combined_df_2013[combined_df_2013$Whale == "F13066",],aes(x= rev(sample_date), y= dN), size = point_size, shape = 21, fill = "#FCBBA1") + 
  geom_line(data= combined_df_2013[combined_df_2013$Whale == "F13068",], aes(x= rev(sample_date), y= dN),col="#FC9272") +
  geom_point(data= combined_df_2013[combined_df_2013$Whale == "F13068",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#FC9272") +
  geom_line(data= combined_df_2013[combined_df_2013$Whale == "F13073",], aes(x= rev(sample_date), y= dN),col="#FB6A4A")+ 
  geom_point(data= combined_df_2013[combined_df_2013$Whale == "F13073",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#FB6A4A") + 
  geom_line(data= combined_df_2013[combined_df_2013$Whale == "F13076",], aes(x= rev(sample_date), y= dN),col="#DE2D26")+ 
  geom_point(data= combined_df_2013[combined_df_2013$Whale == "F13076",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#DE2D26") +
  geom_line(data= combined_df_2013[combined_df_2013$Whale == "F13083",], aes(x= rev(sample_date), y= dN),col="#A50F15")+ 
  geom_point(data= combined_df_2013[combined_df_2013$Whale == "F13083",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#A50F15") +
  theme_article(base_size = 10) + 
  xlab("Date") +
  ylab(expression(paste(delta^{15}, "N (\u2030)")))  #+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#Relació de la senyal isotòpica de carboni amb la data estimada de cada punt
ggplot() + geom_line(data= combined_df_2013[combined_df_2013$Whale == "F13065",], aes(x= rev(sample_date), y= dC),col="#FEE5D9")+ 
  geom_point(data= combined_df_2013[combined_df_2013$Whale == "F13065",],aes(x= rev(sample_date), y= dC), size = point_size, shape = 21, fill = "#FEE5D9") +
  geom_line(data= combined_df_2013[combined_df_2013$Whale == "F13066",], aes(x= rev(sample_date), y= dC),col="#FCBBA1")+ 
  geom_point(data= combined_df_2013[combined_df_2013$Whale == "F13066",],aes(x= rev(sample_date), y= dC), size = point_size, shape = 21, fill = "#FCBBA1") + 
  geom_line(data= combined_df_2013[combined_df_2013$Whale == "F13068",], aes(x= rev(sample_date), y= dC),col="#FC9272") +
  geom_point(data= combined_df_2013[combined_df_2013$Whale == "F13068",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#FC9272") +
  geom_line(data= combined_df_2013[combined_df_2013$Whale == "F13073",], aes(x= rev(sample_date), y= dC),col="#FB6A4A")+ 
  geom_point(data= combined_df_2013[combined_df_2013$Whale == "F13073",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#FB6A4A") + 
  geom_line(data= combined_df_2013[combined_df_2013$Whale == "F13076",], aes(x= rev(sample_date), y= dC),col="#DE2D26")+ 
  geom_point(data= combined_df_2013[combined_df_2013$Whale == "F13076",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#DE2D26") +
  geom_line(data= combined_df_2013[combined_df_2013$Whale == "F13083",], aes(x= rev(sample_date), y= dC),col="#A50F15")+ 
  geom_point(data= combined_df_2013[combined_df_2013$Whale == "F13083",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#A50F15") +
  theme_article(base_size = 10) + 
  xlab("Date") +
  ylab(expression(paste(delta^{13}, "C (\u2030)")))  #+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#MEAN 2013
merge_2013_filtered <- combined_df_2013 %>% dplyr::group_by(Cm) %>% dplyr::summarise(across(.cols = where(is.numeric), .fns = ~mean(dN, na.rm = TRUE))) #calcular la mitjana per centímetre ja que suposem que estàn alineats
ggplot(merge_2013_filtered, aes(x=merge_2013_filtered$Cm,y=merge_2013_filtered$dN)) + geom_point() + geom_line() + theme_article() + xlab("Cm") + ylab(expression(paste(delta^{15}, "N (\u2030)")))  
merge_2013_filtered <- combined_df_2013 %>% dplyr::group_by(Cm) %>% dplyr::summarise(across(.cols = where(is.numeric), .fns = ~mean(dC, na.rm = TRUE))) #calcular la mitjana per centímetre ja que suposem que estàn alineats
ggplot(merge_2013_filtered, aes(x=merge_2013_filtered$Cm,y=merge_2013_filtered$dC)) + geom_point() + geom_line() + theme_article() + xlab("Cm") + ylab(expression(paste(delta^{13}, "C (\u2030)")))  

#________________________________________________________________________________
#Mostres de l'any 2015 (n=9)
#________________________________________________________________________________

whale_isos_2015 <- dplyr::filter(whale_isos, Year == 2015)
unique(whale_isos_2015$Whale) #(n=9)

point_size <- 2
combined_df_2015 <- NULL;
for (i in 1:length(unique(whale_isos_2015$Whale))) {
  df_i <- dplyr::filter(whale_isos, Whale == unique(whale_isos_2015$Whale)[i]) #Filtrem l'individu que volem
  growth_rate <- 16  # baleen growth rate in centimetres per year
  last_sample <- unique(ymd(df_i$Data_capt)) # final date of last keratin sample
  df_i <- 
    df_i %>%
    mutate(days = Cm * 365.25 / growth_rate) %>%
    mutate(days = days - days[1]) %>%
    mutate(rev_days = days - days[length(days)]) %>%
    mutate(sample_date = last_sample + rev_days) %>%
    mutate(year = sample_date)%>%
    mutate(year_rev = rev(sample_date)) #Funció que retrospectivament data cada punt analitzat
  
  df_i$Year_from_sample_date <- format(as.Date(df_i$year_rev), "%Y") #extreure l'any de la data estimada
  df_i$Year_month <- format(as.Date(df_i$year_rev), "%Y/%m") #extreure l'any i el mes de la data estimada
  df_i$Month_from_sample_date <- format(as.Date(df_i$year_rev), "%m") #extreure només el mes de la data estimada
  isop <- ggplot(df_i, aes(x = rev(Cm))) + 
    theme_classic(base_size = 16) + 
    xlab("Baleen sample (cm from youngest sample)") + 
    scale_x_continuous(breaks = df_i$Cm, 
                       labels = rev(df_i$Cm),
                       sec.axis = sec_axis(~., 
                                           breaks = rev(df_i$Cm),
                                           labels = rev(df_i$sample_date),
                                           name = "Date"))
  #Aquest pas es pot saltar si es vol obtenir un plot amb doble y-axis
  isop <- ggplot(df_i, aes(x = rev(sample_date))) + 
    xlab("Date") + 
    scale_x_continuous(breaks = seq(min(df_i$sample_date), max(df_i$sample_date), by = "1 month"),
                       labels = format(seq(min(rev(df_i$sample_date)), max(rev(df_i$sample_date)), by = "1 month"), "%b-%Y"))
  isop <- isop + geom_line(aes(y = dC),col = "black") + 
    geom_point(aes(y = dC), size = point_size, shape = 21, fill = "black", col = "black") +
    ylab(expression(paste(delta^{13}, "C (\u2030)")))+ theme_article()  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  a <- mean(na.omit(df_i$dN))# add the d15N data but need to scale it to d13C
  a
  b <- mean(na.omit(df_i$dC))# add the d15N data but need to scale it to d13C
  b
  v_shift <- a - b
  isop <- isop + geom_line(aes(y = dN - v_shift), col = "darkgrey") + 
    geom_point(aes(y = (dN - v_shift)), size = point_size, 
               shape = 21, fill = "darkgrey", col = "darkgrey")
  isop <- isop + 
    scale_y_continuous(sec.axis = sec_axis(~.+v_shift, 
                                           name = expression(paste(delta^{15}, "N (\u2030)"))),limits = c(-30, -10)) + labs(title= unique(whale_isos_2015$Whale)[i])
  new_df <- df_i
  combined_df_2015 <- rbind(combined_df_2015, new_df)
  print(isop)
}

library(RColorBrewer)
brewer.pal(n = 9, name = "Greens")
#Relació de la senyal isotòpica de nitrogen amb la data estimada de cada punt
ggplot() + geom_line(data= combined_df_2015[combined_df_2015$Whale == "F15078",], aes(x= rev(sample_date), y= dN),col="#F7FCF5")+ 
  geom_point(data= combined_df_2015[combined_df_2015$Whale == "F15078",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#F7FCF5") + 
  geom_line(data= combined_df_2015[combined_df_2015$Whale == "F15097",], aes(x= rev(sample_date), y= dN),col="#E5F5E0") +
  geom_point(data= combined_df_2015[combined_df_2015$Whale == "F15097",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#E5F5E0") +
  geom_line(data= combined_df_2015[combined_df_2015$Whale == "F15086",], aes(x= rev(sample_date), y= dN),col="#C7E9C0")+ 
  geom_point(data= combined_df_2015[combined_df_2015$Whale == "F15086",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#C7E9C0") + 
  geom_line(data= combined_df_2015[combined_df_2015$Whale == "F15083",], aes(x= rev(sample_date), y= dN),col="#A1D99B")+ 
  geom_point(data= combined_df_2015[combined_df_2015$Whale == "F15083",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#A1D99B") +
  geom_line(data= combined_df_2015[combined_df_2015$Whale == "F15080",], aes(x= rev(sample_date), y= dN),col="#74C476")+ 
  geom_point(data= combined_df_2015[combined_df_2015$Whale == "F15080",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#74C476") + 
  geom_line(data= combined_df_2015[combined_df_2015$Whale == "F15079",], aes(x= rev(sample_date), y= dN),col="#41AB5D") +
  geom_point(data= combined_df_2015[combined_df_2015$Whale == "F15079",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#41AB5D") +
  geom_line(data= combined_df_2015[combined_df_2015$Whale == "F15088",], aes(x= rev(sample_date), y= dN),col="#238B45")+ 
  geom_point(data= combined_df_2015[combined_df_2015$Whale == "F15088",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#238B45") +
  geom_line(data= combined_df_2015[combined_df_2015$Whale == "F15084",], aes(x= rev(sample_date), y= dN),col="#006D2C")+ 
  geom_point(data= combined_df_2015[combined_df_2015$Whale == "F15084",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#006D2C") + 
  geom_line(data= combined_df_2015[combined_df_2015$Whale == "F15087",], aes(x= rev(sample_date), y= dN),col="#00441B")+ 
  geom_point(data= combined_df_2015[combined_df_2015$Whale == "F15087",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#00441B") +
  theme_article(base_size = 10) + 
  xlab("Date") +
  ylab(expression(paste(delta^{15}, "N (\u2030)")))  #+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#Relació de la senyal isotòpica de carboni amb la data estimada de cada punt
ggplot() + geom_line(data= combined_df_2015[combined_df_2015$Whale == "F15078",], aes(x= rev(sample_date), y= dC),col="#F7FCF5")+ 
  geom_point(data= combined_df_2015[combined_df_2015$Whale == "F15078",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#F7FCF5") + 
  geom_line(data= combined_df_2015[combined_df_2015$Whale == "F15097",], aes(x= rev(sample_date), y= dC),col="#E5F5E0") +
  geom_point(data= combined_df_2015[combined_df_2015$Whale == "F15097",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#E5F5E0") +
  geom_line(data= combined_df_2015[combined_df_2015$Whale == "F15086",], aes(x= rev(sample_date), y= dC),col="#C7E9C0")+ 
  geom_point(data= combined_df_2015[combined_df_2015$Whale == "F15086",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#C7E9C0") + 
  geom_line(data= combined_df_2015[combined_df_2015$Whale == "F15083",], aes(x= rev(sample_date), y= dC),col="#A1D99B")+ 
  geom_point(data= combined_df_2015[combined_df_2015$Whale == "F15083",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#A1D99B") +
  geom_line(data= combined_df_2015[combined_df_2015$Whale == "F15080",], aes(x= rev(sample_date), y= dC),col="#74C476")+ 
  geom_point(data= combined_df_2015[combined_df_2015$Whale == "F15080",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#74C476") + 
  geom_line(data= combined_df_2015[combined_df_2015$Whale == "F15079",], aes(x= rev(sample_date), y= dC),col="#41AB5D") +
  geom_point(data= combined_df_2015[combined_df_2015$Whale == "F15079",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#41AB5D") +
  geom_line(data= combined_df_2015[combined_df_2015$Whale == "F15088",], aes(x= rev(sample_date), y= dC),col="#238B45")+ 
  geom_point(data= combined_df_2015[combined_df_2015$Whale == "F15088",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#238B45") +
  geom_line(data= combined_df_2015[combined_df_2015$Whale == "F15084",], aes(x= rev(sample_date), y= dC),col="#006D2C")+ 
  geom_point(data= combined_df_2015[combined_df_2015$Whale == "F15084",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#006D2C") + 
  geom_line(data= combined_df_2015[combined_df_2015$Whale == "F15087",], aes(x= rev(sample_date), y= dC),col="#00441B")+ 
  geom_point(data= combined_df_2015[combined_df_2015$Whale == "F15087",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#00441B") +
  theme_article(base_size = 10) + 
  xlab("Date") +
  ylab(expression(paste(delta^{13}, "C (\u2030)")))  #+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#MEAN 2015
merge_2015_filtered <- combined_df_2015 %>% dplyr::group_by(Cm) %>% dplyr::summarise(across(.cols = where(is.numeric), .fns = ~mean(dN, na.rm = TRUE))) #calcular la mitjana per centímetre ja que suposem que estàn alineats
ggplot(merge_2015_filtered, aes(x=merge_2015_filtered$Cm,y=merge_2015_filtered$dN)) + geom_point() + geom_line() + theme_article() + xlab("Cm") + ylab(expression(paste(delta^{15}, "N (\u2030)")))  
merge_2015_filtered <- combined_df_2015 %>% dplyr::group_by(Cm) %>% dplyr::summarise(across(.cols = where(is.numeric), .fns = ~mean(dC, na.rm = TRUE))) #calcular la mitjana per centímetre ja que suposem que estàn alineats
ggplot(merge_2015_filtered, aes(x=merge_2015_filtered$Cm,y=merge_2015_filtered$dC)) + geom_point() + geom_line() + theme_article() + xlab("Cm") + ylab(expression(paste(delta^{13}, "C (\u2030)")))  

#________________________________________________________________________________
#Mostres de l'any 2018 (n=12)
#________________________________________________________________________________

whale_isos_2018 <- dplyr::filter(whale_isos, Year == 2018)
unique(whale_isos_2018$Whale) #(n=12)

point_size <- 2
combined_df_2018 <- NULL;
for (i in 1:length(unique(whale_isos_2018$Whale))) {
  df_i <- dplyr::filter(whale_isos, Whale == unique(whale_isos_2018$Whale)[i]) #Filtrem l'individu que volem
  growth_rate <- 16  # baleen growth rate in centimetres per year
  last_sample <- unique(ymd(df_i$Data_capt)) # final date of last keratin sample
  df_i <- 
    df_i %>%
    mutate(days = Cm * 365.25 / growth_rate) %>%
    mutate(days = days - days[1]) %>%
    mutate(rev_days = days - days[length(days)]) %>%
    mutate(sample_date = last_sample + rev_days) %>%
    mutate(year = sample_date)%>%
    mutate(year_rev = rev(sample_date)) #Funció que retrospectivament data cada punt analitzat
  
  df_i$Year_from_sample_date <- format(as.Date(df_i$year_rev), "%Y") #extreure l'any de la data estimada
  df_i$Year_month <- format(as.Date(df_i$year_rev), "%Y/%m") #extreure l'any i el mes de la data estimada
  df_i$Month_from_sample_date <- format(as.Date(df_i$year_rev), "%m") #extreure només el mes de la data estimada
  isop <- ggplot(df_i, aes(x = rev(Cm))) + 
    theme_classic(base_size = 16) + 
    xlab("Baleen sample (cm from youngest sample)") + 
    scale_x_continuous(breaks = df_i$Cm, 
                       labels = rev(df_i$Cm),
                       sec.axis = sec_axis(~., 
                                           breaks = rev(df_i$Cm),
                                           labels = rev(df_i$sample_date),
                                           name = "Date"))
  #Aquest pas es pot saltar si es vol obtenir un plot amb doble y-axis
  isop <- ggplot(df_i, aes(x = rev(sample_date))) + 
    xlab("Date") + 
    scale_x_continuous(breaks = seq(min(df_i$sample_date), max(df_i$sample_date), by = "1 month"),
                       labels = format(seq(min(rev(df_i$sample_date)), max(rev(df_i$sample_date)), by = "1 month"), "%b-%Y"))
  isop <- isop + geom_line(aes(y = dC),col = "black") + 
    geom_point(aes(y = dC), size = point_size, shape = 21, fill = "black", col = "black") +
    ylab(expression(paste(delta^{13}, "C (\u2030)")))+ theme_article()  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  a <- mean(na.omit(df_i$dN))# add the d15N data but need to scale it to d13C
  a
  b <- mean(na.omit(df_i$dC))# add the d15N data but need to scale it to d13C
  b
  v_shift <- a - b
  isop <- isop + geom_line(aes(y = dN - v_shift), col = "darkgrey") + 
    geom_point(aes(y = (dN - v_shift)), size = point_size, 
               shape = 21, fill = "darkgrey", col = "darkgrey")
  isop <- isop + 
    scale_y_continuous(sec.axis = sec_axis(~.+v_shift, 
                                           name = expression(paste(delta^{15}, "N (\u2030)"))),limits = c(-30, -10)) + labs(title= unique(whale_isos_2018$Whale)[i])
  new_df <- df_i
  combined_df_2018 <- rbind(combined_df_2018, new_df)
  print(isop)
}

library(RColorBrewer)
brewer.pal(n = 9, name = "Oranges")
#Relació de la senyal isotòpica de nitrogen amb la data estimada de cada punt
ggplot() + geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18003",], aes(x= rev(sample_date), y= dN),col="#FFF5EB")+ 
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18003",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#FFF5EB") + 
  geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18007",], aes(x= rev(sample_date), y= dN),col="#FEE6CE") +
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18007",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#FEE6CE") +
  geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18009",], aes(x= rev(sample_date), y= dN),col="#FDD0A2")+ 
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18009",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#FDD0A2") + 
  geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18030",], aes(x= rev(sample_date), y= dN),col="#FDAE6B")+ 
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18030",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#FDAE6B") +
  geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18038",], aes(x= rev(sample_date), y= dN),col="#FD8D3C")+ 
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18038",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#FD8D3C") + 
  geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18043",], aes(x= rev(sample_date), y= dN),col="#F16913") +
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18043",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#F16913") +
  geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18044",], aes(x= rev(sample_date), y= dN),col="#D94801")+ 
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18044",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#D94801") +
  geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18050",], aes(x= rev(sample_date), y= dN),col="#A63603")+ 
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18050",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#A63603") +
  geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18051",], aes(x= rev(sample_date), y= dN),col="#7F2704")+ 
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18051",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#7F2704") + 
  geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18054",], aes(x= rev(sample_date), y= dN),col="brown4") +
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18054",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "brown4") +
  #geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18098",], aes(x= rev(sample_date), y= dN),col="#AD4523FF")+ 
  #geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18098",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#AD4523FF") +
  #geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18025",], aes(x= rev(sample_date), y= dN),col="#9E3D22FF")+ 
  #geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18025",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#9E3D22FF") +
  #geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18022",], aes(x= rev(sample_date), y= dN),col="#9E3D22FF")+ 
  #geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18022",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#9E3D22FF") +
  theme_article(base_size = 10) + 
  xlab("Date") +
  ylab(expression(paste(delta^{15}, "N (\u2030)")))  #+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#Relació de la senyal isotòpica de carboni amb la data estimada de cada punt
ggplot() + geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18003",], aes(x= rev(sample_date), y= dC),col="#FFF5EB")+ 
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18003",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#FFF5EB") + 
  geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18007",], aes(x= rev(sample_date), y= dC),col="#FEE6CE") +
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18007",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#FEE6CE") +
  geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18009",], aes(x= rev(sample_date), y= dC),col="#FDD0A2")+ 
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18009",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#FDD0A2") + 
  geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18030",], aes(x= rev(sample_date), y= dC),col="#FDAE6B")+ 
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18030",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#FDAE6B") +
  geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18038",], aes(x= rev(sample_date), y= dC),col="#FD8D3C")+ 
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18038",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#FD8D3C") + 
  geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18043",], aes(x= rev(sample_date), y= dC),col="#F16913") +
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18043",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#F16913") +
  geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18044",], aes(x= rev(sample_date), y= dC),col="#D94801")+ 
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18044",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#D94801") +
  geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18050",], aes(x= rev(sample_date), y= dC),col="#A63603")+ 
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18050",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#A63603") +
  geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18051",], aes(x= rev(sample_date), y= dC),col="#7F2704")+ 
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18051",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#7F2704") + 
  geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18054",], aes(x= rev(sample_date), y= dC),col="brown4") +
  geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18054",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "brown4") +
  #geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18098",], aes(x= rev(sample_date), y= dC),col="#AD4523FF")+ 
  #geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18098",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#AD4523FF") +
  #geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18025",], aes(x= rev(sample_date), y= dC),col="#9E3D22FF")+ 
  #geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18025",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#9E3D22FF") +
  #geom_line(data= combined_df_2018[combined_df_2018$Whale == "F18022",], aes(x= rev(sample_date), y= dC),col="#9E3D22FF")+ 
  #geom_point(data= combined_df_2018[combined_df_2018$Whale == "F18022",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#9E3D22FF") +
  theme_article(base_size = 10) + 
  xlab("Date") +
  ylab(expression(paste(delta^{13}, "C (\u2030)")))  #+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#MEAN 2018
merge_2018_filtered <- combined_df_2018 %>% dplyr::group_by(Cm) %>% dplyr::summarise(across(.cols = where(is.numeric), .fns = ~mean(dN, na.rm = TRUE))) #calcular la mitjana per centímetre ja que suposem que estàn alineats
ggplot(merge_2018_filtered, aes(x=merge_2018_filtered$Cm,y=merge_2018_filtered$dN)) + geom_point() + geom_line() + theme_article() + xlab("Cm") + ylab(expression(paste(delta^{15}, "N (\u2030)")))  

merge_2018_filtered <- combined_df_2018 %>% dplyr::group_by(Cm) %>% dplyr::summarise(across(.cols = where(is.numeric), .fns = ~mean(dC, na.rm = TRUE))) #calcular la mitjana per centímetre ja que suposem que estàn alineats
ggplot(merge_2018_filtered, aes(x=merge_2018_filtered$Cm,y=merge_2018_filtered$dC)) + geom_point() + geom_line() + theme_article() + xlab("Cm") + ylab(expression(paste(delta^{13}, "C (\u2030)")))  

#________________________________________________________________________________
#Mostres de l'any 2022 (n=4)
#________________________________________________________________________________

whale_isos_2022 <- dplyr::filter(whale_isos, Year == 2022)
unique(whale_isos_2022$Whale) #(n=4)

require(scales)
point_size <- 3
combined_df_2022 <- NULL;
for (i in 1:length(unique(whale_isos_2022$Whale))) {
  df_i <- dplyr::filter(whale_isos, Whale == unique(whale_isos_2022$Whale)[i]) #Filtrem l'individu que volem
  growth_rate <- 16  # baleen growth rate in centimetres per year
  last_sample <- unique(ymd(df_i$Data_capt)) # final date of last keratin sample
  
  
  df_i <- 
    df_i %>%
    mutate(days = Cm * 365.25 / growth_rate) %>%
    mutate(days = days - days[1]) %>%
    mutate(rev_days = days - days[length(days)]) %>%
    mutate(sample_date = last_sample + rev_days) %>%
    mutate(year = sample_date)%>%
    mutate(year_rev = rev(sample_date)) #Funció que retrospectivament data cada punt analitzat
  
  df_i$Year_from_sample_date <- format(as.Date(df_i$year_rev), "%Y") #extreure l'any de la data estimada
  df_i$Year_month <- format(as.Date(df_i$year_rev), "%Y/%m") #extreure l'any i el mes de la data estimada
  df_i$Month_from_sample_date <- format(as.Date(df_i$year_rev), "%m") #extreure només el mes de la data estimada
  isop <- ggplot(df_i, aes(x = rev(Cm))) + 
    theme_bw(base_size = 20) + 
    xlab("Baleen sample (cm from youngest sample)") + 
    scale_x_continuous(breaks = df_i$Cm, 
                       labels = rev(df_i$Cm),
                       sec.axis = sec_axis(~., 
                                           breaks = rev(df_i$Cm),
                                           labels = rev(df_i$sample_date),
                                           name = "Date"))
  #Aquest pas es pot saltar si es vol obtenir un plot amb doble y-axis
  isop <- ggplot(df_i, aes(x = rev(sample_date))) + 
    xlab("Date") + scale_x_date(labels = scales::date_format("%Y"))
  
  # Add δ13C with a color aesthetic
  isop <- isop + 
    geom_line(aes(y = dC, color = "delta13C"), col = "black") + 
    geom_point(aes(y = dC, color = "delta13C"), size = point_size, shape = 21, fill = "black", col = "black") +
    ylab(expression(paste(delta^{13}, "C (\u2030)")))
  
  # Add δ15N with a color aesthetic and vertical shift
  isop <- isop + 
    geom_line(aes(y = dN - v_shift, color = "delta15N"), col = "darkgrey") + 
    geom_point(aes(y = (dN - v_shift), color = "delta15N"), size = point_size, shape = 21, fill = "darkgrey", col = "darkgrey")
  
  # Add the legend for the stable isotopes
  isop <- isop + 
    scale_color_manual(values = c("delta13C" = "black", "delta15N" = "darkgrey"),
                       name = "Stable Isotopes",
                       labels = c("δ13C", "δ15N")) +
    theme_bw(base_size = 20)
  
  
  a <- mean(na.omit(df_i$dN))# add the d15N data but need to scale it to d13C
  a
  b <- mean(na.omit(df_i$dC))# add the d15N data but need to scale it to d13C
  b
  v_shift <- a - b
  isop <- isop + geom_line(aes(y = dN - v_shift), col = "darkgrey") + 
    geom_point(aes(y = (dN - v_shift)), size = point_size, 
               shape = 21, fill = "darkgrey", col = "darkgrey")
  isop <- isop + 
    scale_y_continuous(sec.axis = sec_axis(~.+v_shift, 
                                           name = expression(paste(delta^{15}, "N (\u2030)"))),limits = c(-22, -17)) #+ labs(title= unique(whale_isos_2022$Whale)[i])
  
  
  new_df <- df_i
  combined_df_2022 <- rbind(combined_df_2022, new_df)
  print(isop)
}


#Relació de la senyal isotòpica de nitrogen amb la data estimada de cada punt
ggplot() + geom_line(data= combined_df_2022[combined_df_2022$Whale == "F22046",], aes(x= rev(sample_date), y= dN),col="#B9DDF1FF")+ 
  geom_point(data= combined_df_2022[combined_df_2022$Whale == "F22046",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#B9DDF1FF") + 
  geom_line(data= combined_df_2022[combined_df_2022$Whale == "F22053",], aes(x= rev(sample_date), y= dN),col="#7EAED3FF") +
  geom_point(data= combined_df_2022[combined_df_2022$Whale == "F22053",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#7EAED3FF") +
  geom_line(data= combined_df_2022[combined_df_2022$Whale == "F22054",], aes(x= rev(sample_date), y= dN),col="#5081AEFF")+ 
  geom_point(data= combined_df_2022[combined_df_2022$Whale == "F22054",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#5081AEFF") + 
  geom_line(data= combined_df_2022[combined_df_2022$Whale == "F22056",], aes(x= rev(sample_date), y= dN),col="#2A5783FF")+ 
  geom_point(data= combined_df_2022[combined_df_2022$Whale == "F22056",],aes(x= rev(sample_date), y = dN), size = point_size, shape = 21, fill = "#2A5783FF") +
  theme_article(base_size = 10) + 
  xlab("Date") +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) #+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#Relació de la senyal isotòpica de carboni amb la data estimada de cada punt
ggplot() + geom_line(data= combined_df_2022[combined_df_2022$Whale == "F22046",], aes(x= rev(sample_date), y= dC),col="#B9DDF1FF")+ 
  geom_point(data= combined_df_2022[combined_df_2022$Whale == "F22046",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#B9DDF1FF") + 
  geom_line(data= combined_df_2022[combined_df_2022$Whale == "F22053",], aes(x= rev(sample_date), y= dC),col="#7EAED3FF") +
  geom_point(data= combined_df_2022[combined_df_2022$Whale == "F22053",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#7EAED3FF") +
  geom_line(data= combined_df_2022[combined_df_2022$Whale == "F22054",], aes(x= rev(sample_date), y= dC),col="#5081AEFF")+ 
  geom_point(data= combined_df_2022[combined_df_2022$Whale == "F22054",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#5081AEFF") + 
  #geom_line(data= combined_df_2022[combined_df_2022$Whale == "F22056",], aes(x= rev(sample_date), y= dC),col="#2A5783FF")+ 
  #geom_point(data= combined_df_2022[combined_df_2022$Whale == "F22056",],aes(x= rev(sample_date), y = dC), size = point_size, shape = 21, fill = "#2A5783FF") +
  theme_article(base_size = 10) + 
  xlab("Date") +
  ylab(expression(paste(delta^{13}, "C (\u2030)"))) #+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#MEAN 2022
merge_2022_filtered <- combined_df_2022 %>% dplyr::group_by(Cm) %>% dplyr::summarise(across(.cols = where(is.numeric), .fns = ~mean(dN, na.rm = TRUE))) #calcular la mitjana per centímetre ja que suposem que estàn alineats
ggplot(merge_2022_filtered, aes(x=merge_2022_filtered$Cm,y=merge_2022_filtered$dN)) + geom_point() + geom_line() + theme_article() + xlab("Cm") + ylab(expression(paste(delta^{15}, "N (\u2030)")))  

merge_2022_filtered <- combined_df_2022 %>% dplyr::group_by(Cm) %>% dplyr::summarise(across(.cols = where(is.numeric), .fns = ~mean(dC, na.rm = TRUE))) #calcular la mitjana per centímetre ja que suposem que estàn alineats
ggplot(merge_2022_filtered, aes(x=merge_2022_filtered$Cm,y=merge_2022_filtered$dC)) + geom_point() + geom_line() + theme_article() + xlab("Cm") + ylab(expression(paste(delta^{13}, "C (\u2030)")))  
#__________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

########################################################
#Merge the individuals from 2013,2015,2018 and 2022 data
########################################################

merge <- rbind(combined_df_2013,combined_df_2015,combined_df_2018, combined_df_2022)

merge <-merge[complete.cases(merge$dC), ] #remove the rows with no dC values in order to be able to correct for the Suess effect


ggplot(merge, aes(as.Date(year_rev),dN, color=Whale)) + geom_line() + theme_classic() 
# Assuming your data frame is named 'merge'
library(ggplot2)

data <- read_excel("~/Desktop/Paper_climate/All_merged_2013_to_2022_Suess_20_oct_df.xlsx")
data$Class = factor(data$Class, levels=c("dN","dC")) #labels = c(expression(paste(delta^{15}, "N (\u2030)")), expression(paste(delta^{13}, "C (\u2030)"))

nb.cols <- 29
mycolors <- colorRampPalette(brewer.pal(8, "RdYlBu"))(nb.cols) #"YlGnBu"

cyl_names <- c(
  "dN" = "delta^{15}*N",
  "dC" = "delta^{13}*C"
)

#stable isotope ratios
rect <- data.frame(xmin=as.Date("2018-07-01"), xmax=as.Date("2019-06-01"), ymin=-Inf, ymax=Inf)
ggplot(data=data, aes(x = as.Date(year_rev), y= Isotope, color = Whale)) +
  geom_smooth(method = "loess", se = FALSE, span = 0.3, aes(fill = Whale))  + facet_wrap(vars(Class), ncol=1, scales = "free", labeller= labeller(Class = as_labeller(cyl_names, label_parsed))) +  # Fill color for geom_rect
  xlab("Year") +
  ylab(expression(paste("Isotope ratio (\u2030)"))) +
  scale_x_date(limits = c(min(merge$year_rev), max(merge$year_rev)+365), date_breaks = "2 year", date_labels = "%Y") +
  geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),fill="white",alpha=1,inherit.aes = FALSE) +
  geom_vline(xintercept = rect$xmin, linetype="dashed", 
             color = "darkgrey") +
  geom_vline(xintercept = rect$xmax, linetype="dashed", 
             color = "darkgrey") + labs(fill= "Whale ID", color= "Whale ID") +
  scale_fill_manual(values = mycolors) +
  scale_color_manual(values = mycolors) + theme(aspect.ratio = 1.5/4, legend.key.size = unit(0.5, 'cm'), #change legend key size
                                                                                legend.key.height = unit(0.5, 'cm'), #change legend key height
                                                                                legend.key.width = unit(0.5, 'cm'), #change legend key width                                                                          legend.title = element_text(size=10), #change legend title font size
                                                                               legend.text = element_text(size=12))




#QUIK PLOT
library(egg)
ggplot(data, aes(as.factor(Year_from_sample_date),Isotope)) + 
  geom_boxplot(fill="#74A9CF", alpha=0.9, width=0.2) +
  geom_jitter(color="black", size=0.1, alpha=0.1, shape=1, width = 0.1) +
  facet_wrap(vars(Class), ncol=1, scales = "free", labeller= labeller(Class = as_labeller(cyl_names, label_parsed))) +  # Fill color for geom_rect
  #stat_summary(fun.y=mean, geom="point", shape=20, size=2, color="red", fill="red") +
  xlab("Year") + 
  ylab(expression(paste("Isotope ratio (\u2030)"))) +
  theme_article(base_size = 10) + theme(aspect.ratio=2/4)      

# Compute the analysis of variance
final_fw.data_d13cor$Year_from_sample_date <- as.factor(final_fw.data_d13cor$Year_from_sample_date)
res.aov <- aov(d13c.cor ~ Year_from_sample_date, data = final_fw.data_d13cor)
# Summary of the analysis
summary(res.aov)

TukeyHSD(res.aov)

nb.cols <- 29
mycolors <- colorRampPalette(brewer.pal(8, "RdYlBu"))(nb.cols) #"YlGnBu"

rect <- data.frame(xmin=as.Date("2018-07-01"), xmax=as.Date("2019-06-01"), ymin=-Inf, ymax=Inf)
ggplot(merge2, aes(x = as.Date(year_rev), y=dN, color = Whale)) +
  geom_smooth(method = "loess", se = FALSE, span = 0.3, aes(fill = Whale)) +  # Fill color for geom_rect
  theme_article(base_size = 20) +
  xlab("Year") +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  scale_x_date(limits = c(min(merge$year_rev), max(merge$year_rev)+365), date_breaks = "2 year", date_labels = "%Y") +
  geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),fill="white",alpha=1,inherit.aes = FALSE) +
  geom_vline(xintercept = rect$xmin, linetype="dashed", 
             color = "darkgrey") +
  geom_vline(xintercept = rect$xmax, linetype="dashed", 
             color = "darkgrey") +
  scale_fill_manual(values = mycolors) +
  scale_color_manual(values = mycolors) + theme(aspect.ratio = 3/4)



rect <- data.frame(xmin=as.Date("2018-07-01"), xmax=as.Date("2019-06-01"), ymin=-Inf, ymax=Inf)
ggplot(merge2, aes(x = as.Date(year_rev), y = dC, color = Whale)) +
  geom_smooth(method = "loess", se = FALSE, span=0.3) +  # Add LOESS line without confidence interval
  theme_article(base_size = 20) +
  xlab("Year") + 
  ylab(expression(paste(delta^{13}, "C (\u2030)"))) +  
  scale_x_date(limits = c(min(merge$year_rev), max(merge$year_rev)+365), date_breaks = "2 year", date_labels = "%Y")+ 
  geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),fill="white",alpha=1,inherit.aes = FALSE) +
  geom_vline(xintercept = rect$xmin, linetype="dashed", 
             color = "darkgrey") +
  geom_vline(xintercept = rect$xmax, linetype="dashed", 
             color = "darkgrey") +
  scale_fill_manual(values = mycolors) +
  scale_color_manual(values = mycolors) + theme(aspect.ratio = 3/4)

rect <- data.frame(xmin=as.Date("2018-07-01"), xmax=as.Date("2019-06-01"), ymin=-Inf, ymax=Inf)
ggplot(fw.data_d13cor, aes(x = as.Date(year_rev), y = d13c.cor, color = Whale)) +
  geom_smooth(method = "loess", se = FALSE, span=0.3) +  # Add LOESS line without confidence interval
  theme_article(base_size = 20) +
  xlab("Year") + 
  ylab(expression(paste(delta^{13}, "C (\u2030)"))) +  
  scale_x_date(limits = c(min(merge$year_rev), max(merge$year_rev)+365), date_breaks = "2 year", date_labels = "%Y")+ 
  geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),fill="white",alpha=1,inherit.aes = FALSE) +
  geom_vline(xintercept = rect$xmin, linetype="dashed", 
             color = "darkgrey") +
  geom_vline(xintercept = rect$xmax, linetype="dashed", 
             color = "darkgrey") +
  scale_fill_manual(values = mycolors) +
  scale_color_manual(values = mycolors) + theme(aspect.ratio = 3/4)


#__________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
############################
#Individual removal from the analysis
############################
#merge1 <- filter(merge, Whale %in%  c("F13065", "F13066", "F13068", "F13073", "F13076", "F13129", 
                                      "F13083", "F15078", "F15097", "F15086", "F15083", "F15080", "F15079", 
                                      "F15084", "F15087", "F15088", "F18044", "F18043", "F18038", "F18050", 
                                      "F18051", "F18030", "F18054", "F18007", "F18003", "F18009", "F18098", 
                                      "F18022", "F22046", "F22053", "F22054", "F22056"))

merge1 <- merge
#Remove hybrids
library(dplyr)
merge1 <- merge1 %>% dplyr::mutate(Species = ifelse(Whale %in% c("F13129", "F18022", "F18098"), "hybrid", "fin")) #trec hibrids
#merge2 <- dplyr::filter(merge1, Species == "fin")
length(unique(merge2$Whale))

#Remove immatures (erratic migration movements) NO CAL
#merge2 <- dplyr::filter(merge2, Status == "immatures")
#__________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

###
#Calcul del efecte de Suess 2022
###
#merge_2022 <- dplyr::filter(merge2, Year_from_sample_date == "2022") #creo un dataset de l'any 2022 per tal de poder corregir manualment l'efecte de Suess

#merge_2022 <- data.frame(id = seq(1385,1421,1), merge_2022)

#merge_2022 <- merge_2022 %>%
  add_column(region = "Subpolar North Atlantic")

#Rename fw.data to df1
#df1 <- merge_2022
#Select the fw.data column names (for RSuess) and name it subset
#subset <- merge_2022[c("dC","Year_from_sample_date","region")]
#Renamee subset columnames
#names(subset) <- c("d13c","year","region")
#Determine year as a numeric variable
#subset$year <- as.numeric(subset$year)
#Determine subset as a dataframe
#subset <- as.data.frame(subset)
#subset <- subset %>% add_column(year.y = 2022)
#subset$d13c.uncor <- subset$d13c
#subset <- subset %>% add_column(Laws.cor = 0.043)
#subset <- subset %>% add_column(Suess.cor = 1.302)
#subset <- subset %>% add_column(net.cor = 1.346)
#subset <- subset %>% dplyr::mutate(d13c.cor = dplyr::select(., c(d13c.uncor,net.cor)) %>% rowSums(na.rm = TRUE))
#df2 <- subset
#merge df1 and df2
#fw.data_d13cor_1 <- cbind(df1,df2)
#fw.data_d13cor_1<-fw.data_d13cor_1[,-c(27:29)]
##
#library(openxlsx) #Guardar
#write.xlsx(merge, file = "All_merged_2013_to_2022.xlsx",sheetName = "dades", append = FALSE)
############################
#Calculate the mean baleen plates isotopic values for all indiviuals (Cal corregir d'alguna manera els labels i breaks dels tres darrers anys)
############################

merge_all_filtered <- merge2 %>% dplyr::group_by(Year_month) %>% dplyr::summarise(across(.cols = where(is.numeric), .fns = ~mean(dN, na.rm = TRUE))) #calcular la mitjana per centímetre ja que suposem que estàn alineats
ggplot() +
  #geom_point(data=merge, aes(x=Year_month,y=dN, group=1)) + 
  geom_line(data=merge2, aes(x=Year_month,y=dN, group=Whale), col="grey") +
    geom_point(data=merge_all_filtered, aes(x=Year_month,y=dN, group=1)) + 
    geom_line(data=merge_all_filtered, aes(x=Year_month,y=dN, group=1)) + 
  theme_classic(base_size=10) + 
  xlab("Time") + 
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) + 
  scale_x_discrete(labels = merge_all_filtered$Year_month[c(6, 18, 30, 42, 54, 66, 78, 90, 105, 117, 129)], breaks = merge_all_filtered$Year_month[c(6, 18, 30, 42, 54, 66, 78, 90, 102, 114, 126)])


merge_all_filtered <- merge2 %>% dplyr::group_by(Year_month) %>% dplyr::summarise(across(.cols = where(is.numeric), .fns = ~mean(dC, na.rm = TRUE))) #calcular la mitjana per centímetre ja que suposem que estàn alineats
ggplot() +
  #geom_point(data=merge, aes(x=Year_month,y=dN, group=1)) + 
  geom_line(data=merge2, aes(x=Year_month,y=dC, group=Whale), col="grey") +
  geom_point(data=merge_all_filtered, aes(x=Year_month,y=dC, group=1)) + 
  geom_line(data=merge_all_filtered, aes(x=Year_month,y=dC, group=1)) + 
  theme_classic(base_size=10) + 
  xlab("Cm") + 
  ylab(expression(paste(delta^{13}, "C (\u2030)"))) + 
  scale_x_discrete(labels = merge_all_filtered$Year_month[c(6, 18, 30, 42, 54, 66, 78, 90, 105, 117, 129)], breaks = merge_all_filtered$Year_month[c(6, 18, 30, 42, 54, 66, 78, 90, 102, 114, 126)])

#-------------------------------------------------------------------------------
#Suess effect correction
#-------------------------------------------------------------------------------
#install.packages("SuessR")
library(SuessR)
#Add a sequence number to each row
#merge1 <- merge2
merge1$id <-seq.int(nrow(merge1))
#Add a region to each row
merge1 <- merge1 %>% add_column(region = "Subpolar North Atlantic")
#Rename fw.data to df1
df1 <- merge1
#Select the fw.data column names (for RSuess) and name it subset
subset <- merge1[c("id", "dC","Year_from_sample_date","region")]
#Renamee subset columnames
names(subset) <- c("id", "d13c","year","region")
#Determine year as a numeric variable
subset$year <- as.numeric(subset$year)
#Determine subset as a dataframe
subset <- as.data.frame(subset)
#Correct the Suess effect to the year 2000
df2 <- SuessR(data=subset)
#merge df1 and df2
final_fw.data_d13cor <- merge(df1,df2,by="id") #2022 non corrected



#rename
#names(fw.data_d13cor_2) <- c("id","Cm", "dN" ,"dC","Whale_1","Whale", "Sex","Status","Talla_fetus_.cm.","Sexe_fetus", "Length","Talla_.peus.", "Data_capt","Lat","Long","Edat","Year", "days","rev_days","sample_date", "year","year_rev","Year_from_sample_date","Year_month","Month_from_sample_date","Species","region","year.y","d13c.uncor","Laws.cor","Suess.cor","net.cor","d13c.cor")



#final_fw.data_d13cor <- rbind(fw.data_d13cor_1, fw.data_d13cor_2)
fw.data_d13cor <- final_fw.data_d13cor
#C corrected
fw.data_d13cor$Year_month_date <- zoo::as.Date(zoo::as.yearmon(fw.data_d13cor$Year_month, "%Y/%m"))
fw.data_d13cor$Year_month_date <- as.Date(fw.data_d13cor$Year_month_date, "%Y-%m-%d")
merge_all_filtered <- fw.data_d13cor %>% dplyr::group_by(Year_month_date) %>% dplyr::summarise(across(.cols = where(is.numeric), .fns = ~mean(d13c.cor, na.rm = TRUE))) #calcular la mitjana per centímetre ja que suposem que estàn alineats

rect <- data.frame(xmin=as.Date("2018-07-01"), xmax=as.Date("2019-06-01"), ymin=-Inf, ymax=Inf)

p <- ggplot() +
  #geom_point(data=merge, aes(x=Year_month,y=dN, group=1)) + 
  geom_line(data=fw.data_d13cor, aes(x=Year_month_date,y=d13c.cor, group=Whale), col="darkgrey") +
  geom_point(data=merge_all_filtered, aes(x=Year_month_date,y=d13c.cor, group=1)) + 
  geom_line(data=merge_all_filtered, aes(x=Year_month_date,y=d13c.cor, group=1)) + 
  theme_article() + 
  xlab("Year") + 
  ylab(expression(paste(delta^{13}, "C (\u2030)"))) +  scale_x_date(limits = c(min(fw.data_d13cor$Year_month_date), max(fw.data_d13cor$Year_month_date)), date_breaks = "1 year", date_labels = "%Y") + geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
              fill="lightgrey",
              alpha=1,
              inherit.aes = FALSE)

fw.data_d13cor$Year_month_date <- zoo::as.Date(zoo::as.yearmon(fw.data_d13cor$Year_month, "%Y/%m"))
fw.data_d13cor$Year_month_date <- as.Date(fw.data_d13cor$Year_month_date, "%Y-%m-%d")
merge_all_filtered <- fw.data_d13cor %>% dplyr::group_by(Year_month_date) %>% dplyr::summarise(across(.cols = where(is.numeric), .fns = ~mean(dN, na.rm = TRUE))) #calcular la mitjana per centímetre ja que suposem que estàn alineats
rect <- data.frame(xmin=as.Date("2018-07-01"), xmax=as.Date("2019-06-01"), ymin=-Inf, ymax=Inf)
p1 <- ggplot() +
  #geom_point(data=merge, aes(x=Year_month,y=dN, group=1)) + 
  geom_line(data=fw.data_d13cor, aes(x=Year_month_date,y=dN, group=Whale), col="darkgrey") +
  geom_point(data=merge_all_filtered, aes(x=Year_month_date,y=dN, group=1)) + 
  geom_line(data=merge_all_filtered, aes(x=Year_month_date,y=dN, group=1)) + 
  theme_article() + 
  xlab("Year") + 
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +  
  scale_x_date(limits = c(min(fw.data_d13cor$Year_month_date), max(fw.data_d13cor$Year_month_date)), date_breaks = "1 year", date_labels = "%Y")+ 
  geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),fill="lightgrey",alpha=1,inherit.aes = FALSE)

#install.packages("patchwork")
library(patchwork)
p/p1

library(openxlsx) #Guardar
write.xlsx(fw.data_d13cor, file = "All_merged_2013_to_2022_Suess_5_gen_AMB_HIBRIDS.xlsx")

#-----------------------------------------------------------------
#Alternativa dos al plot de nitrogen i carboni i la mitja
#-----------------------------------------------------------------
#Requires the packages tidyr, ggplot2

require(tidyr)
require(ggplot2)


######################## Set some parameters #########################
#Confidence level, typically 95% = 0.95
Conf_level <-  0.95

#Start and end point of stimulation
stim_start <- "2018-07-25"
stim_end <- "2019-09-25"

######################### import the data ##############################
###
#NITROGEN
###
#Read a text file (comma separated values)
#df_wide <- read.csv("Downloads/FRET-ratio-wide.csv", na.strings = "")

df_wide <- fw.data_d13cor %>%
  select(Whale, dN, year_rev) %>%
  pivot_wider(names_from = Whale, values_from = dN) 

df_wide <- as.data.frame(df_wide)
#Tidy the data, i.e. long format with each row is variable
#df_tidy <- gather(df_wide, Cell, Ratio, -Time)
df_wide$year_rev <- as.Date(df_wide$year_rev)
df_tidy <- gather(df_wide, Whale, dN, -year_rev)
df_tidy$Month_Yr <- format(as.Date(df_tidy$year_rev), "%Y-%m")
library(tidyr)
df_tidy <- df_tidy %>% tidyr::drop_na(dN)
df_tidy <- df_tidy %>% tidyr::drop_na(Whale)
df_tidy <- df_tidy %>% tidyr::drop_na(year_rev)
df_tidy <- as.data.frame(df_tidy)
#Clear the df.summary dataframe if it exist
if (exists("df_summary")) rm(df_summary)


######### Calulcate summary statistics to fill dataframe 'df_summary' ########
# This is tidyverse approach

require(magrittr)
require(dplyr)
library(tidyverse)

df_summary <- df_tidy %>% dplyr::group_by(Month_Yr) %>% dplyr::summarise(mean = mean(dN, na.rm = TRUE),
sd = sd(dN, na.rm = TRUE),
n = n()) 

df_summary <- filter(df_summary, n >1) 


df_summary$Date_half_month <- as.Date(with(df_summary, paste(Month_Yr, "15",sep="-")), "%Y-%m-%d")
######### Calulcate summary statistics to fill dataframe 'df_summary' ########
# This is base R approach

#df_summary <- data.frame(Time=df_wide$year_rev, n=tapply(df_tidy$dN, df_tidy$year_rev, length), mean=tapply(df_tidy$dN, df_tidy$year_rev, mean))

#Add SD and standard error of the mean to the dataframe
#df_summary$sd <- tapply(df_tidy$dN, df_tidy$year_rev, sd)
#df_summary$sem <- df_summary$sd/sqrt(df_summary$n-1)

#Add 95% CI of the mean to the dataframe
#df_summary$CI_lower <- df_summary$mean + qt((1-Conf_level)/2, df=df_summary$n-1)*df_summary$sem
#df_summary$CI_upper <- df_summary$mean - qt((1-Conf_level)/2, df=df_summary$n-1)*df_summary$sem

####################################################
rect <- data.frame(xmin=as.Date("2018-07-01"), xmax=as.Date("2019-06-01"), ymin=-Inf, ymax=Inf)
#### Command to prepare the plot ####
library(dplyr)
library(ggplot2)

# Set your confidence level (e.g., 95%)
Conf_level <- 0.95

# Assuming df_tidy contains your monthly time series data
df_tidy$Mon <- month(df_tidy$year_rev)
df_tidy$Yea <- year(df_tidy$year_rev)
# Group the data by month and year
df_summary <- df_tidy %>%
  dplyr::group_by(Yea, Mon) %>%
  dplyr::summarise(
    mean = mean(dN, na.rm = TRUE),
    CI_lower = quantile(dN, (1 - Conf_level) / 2),
    CI_upper = quantile(dN, 1 - (1 - Conf_level) / 2),
    n = n()
  ) %>%
  filter(n > 1)  # Optional: Remove months with sample size less than 2

# Create a date variable based on Year and Month
df_summary$Date_monthly <- as.Date(paste(df_summary$Yea, df_summary$Mon, "01", sep = "-"))

# Create the ggplot visualization
p1 <- ggplot(df_summary, aes(x = Date_monthly, y = mean)) +
geom_line(color = "#74A9CF", linewidth=2) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), fill = "#74A9CF", alpha = 0.5,color = "black", linetype = "dotted") +
  geom_line(aes(ymin = CI_lower, ymax = CI_upper)) + 
  geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),fill="white",alpha=1,inherit.aes = FALSE) +
  geom_vline(xintercept = rect$xmin, linetype="dashed", 
                                                                color = "darkgrey") +
  geom_vline(xintercept = rect$xmax, linetype="dashed", 
             color = "darkgrey") +
  
 
  # Set x axis labels
  scale_x_date(limits = c(min(df_tidy$year_rev), max(df_tidy$year_rev)), date_breaks = "1 year", date_labels = "%Y")+ 
  
  ## Set the Y-axis scale, remove for autoscale
  
  coord_cartesian(ylim = c(5, 15)) +
  scale_y_continuous(breaks=c(5, 7, 9, 11, 13, 15)) +
  
  ## Set theme&basic font
  theme_article(base_size = 20) +
  
  
  ### Style the axis (font size)
  # theme(axis.text.x = element_text(size=16, angle=0, vjust = 0), axis.text.y = element_text(size=16)) +
  
  ### Set layout of the graph 
  # theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black", fill=NA)) +
  
  ### Set label of x- and y-axis
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) + xlab("Estimated date (years)") +
  
  ### Set aspect ratio of the graph n/n = square
  theme(aspect.ratio=2/4)

p1
  
#p1 <- ggplot(df_summary, aes(x = Date_monthly, y = mean)) +
  geom_line() +
  #geom_line(data=df_tidy, aes(x=as.Date(df_tidy$year_rev), y=dN, group=Whale), color="lightgrey", size=0.25, alpha=0.3) +
  geom_line(size = 0.9, alpha = 0.8, color="#74A9CF") +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), fill = "#74A9CF", color= "#74A9CF", linetype= "dotted", alpha = 0.3) +
  geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),fill="white",alpha=1,inherit.aes = FALSE) +
  geom_vline(xintercept = rect$xmin, linetype="dashed", 
               color = "darkgrey") +
  geom_vline(xintercept = rect$xmax, linetype="dashed", 
             color = "darkgrey") +
  
  
  # Set x axis labels
  scale_x_date(limits = c(min(df_tidy$year_rev), max(df_tidy$year_rev)), date_breaks = "1 year", date_labels = "%Y")+ 
  
  ## Set the Y-axis scale, remove for autoscale
  
  coord_cartesian(ylim = c(5, 15)) +
  scale_y_continuous(breaks=c(5, 7, 9, 11, 13, 15)) +
  
  ## Set theme&basic font
  theme_article(base_size = 20) +
  
  
  ### Style the axis (font size)
  # theme(axis.text.x = element_text(size=16, angle=0, vjust = 0), axis.text.y = element_text(size=16)) +
  
  ### Set layout of the graph 
  # theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black", fill=NA)) +
  
  ### Set label of x- and y-axis
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) + xlab("Estimated date (years)") +
  
  ### Set aspect ratio of the graph n/n = square
  theme(aspect.ratio=2/4)

p1
       

###
#CARBON
###
#Read a text file (comma separated values)
#df_wide <- read.csv("Downloads/FRET-ratio-wide.csv", na.strings = "")

df_wide1 <- fw.data_d13cor %>%
  select(Whale, d13c.cor, year_rev) %>%
  pivot_wider(names_from = Whale, values_from = d13c.cor) 

df_wide1 <- as.data.frame(df_wide1)
#Tidy the data, i.e. long format with each row is variable
#df_tidy <- gather(df_wide, Cell, Ratio, -Time)
df_wide1$year_rev <- as.Date(df_wide1$year_rev)
df_tidy1 <- gather(df_wide1, Whale, d13c.cor, -year_rev)
df_tidy1$Month_Yr <- format(as.Date(df_tidy1$year_rev), "%Y-%m")
library(tidyr)
df_tidy1 <- df_tidy1 %>% tidyr::drop_na(d13c.cor)
df_tidy1 <- df_tidy1 %>% tidyr::drop_na(Whale)
df_tidy1 <- df_tidy1 %>% tidyr::drop_na(year_rev)
df_tidy1 <- as.data.frame(df_tidy1)
#Clear the df.summary dataframe if it exist
if (exists("df_summary")) rm(df_summary)


######### Calulcate summary statistics to fill dataframe 'df_summary' ########
# This is tidyverse approach

require(magrittr)
require(dplyr)
library(tidyverse)

df_summary1 <- df_tidy1 %>% dplyr::group_by(Month_Yr) %>% dplyr::summarise(mean = mean(d13c.cor, na.rm = TRUE),
                                                                         sd = sd(d13c.cor, na.rm = TRUE),
                                                                         n = n()) 

df_summary1 <- filter(df_summary1, n >1) 


df_summary1$Date_half_month <- as.Date(with(df_summary1, paste(Month_Yr, "15",sep="-")), "%Y-%m-%d")
######### Calulcate summary statistics to fill dataframe 'df_summary' ########
# This is base R approach

#df_summary <- data.frame(Time=df_wide$year_rev, n=tapply(df_tidy$dN, df_tidy$year_rev, length), mean=tapply(df_tidy$dN, df_tidy$year_rev, mean))

#Add SD and standard error of the mean to the dataframe
#df_summary$sd <- tapply(df_tidy$dN, df_tidy$year_rev, sd)
#df_summary$sem <- df_summary$sd/sqrt(df_summary$n-1)

#Add 95% CI of the mean to the dataframe
#df_summary$CI_lower <- df_summary$mean + qt((1-Conf_level)/2, df=df_summary$n-1)*df_summary$sem
#df_summary$CI_upper <- df_summary$mean - qt((1-Conf_level)/2, df=df_summary$n-1)*df_summary$sem

####################################################
rect <- data.frame(xmin=as.Date("2018-07-01"), xmax=as.Date("2019-06-01"), ymin=-Inf, ymax=Inf)
#### Command to prepare the plot ####
library(dplyr)
library(ggplot2)

# Set your confidence level (e.g., 95%)
Conf_level <- 0.95

# Assuming df_tidy contains your monthly time series data
df_tidy1$Mon <- month(df_tidy1$year_rev)
df_tidy1$Yea <- year(df_tidy1$year_rev)
# Group the data by month and year
df_summary1 <- df_tidy1 %>%
  dplyr::group_by(Yea, Mon) %>%
  dplyr::summarise(
    mean = mean(d13c.cor, na.rm = TRUE),
    CI_lower = quantile(d13c.cor, (1 - Conf_level) / 2),
    CI_upper = quantile(d13c.cor, 1 - (1 - Conf_level) / 2),
    n = n()
  ) %>%
  filter(n > 1)  # Optional: Remove months with sample size less than 2

# Create a date variable based on Year and Month
df_summary1$Date_monthly <- as.Date(paste(df_summary1$Yea, df_summary1$Mon, "01", sep = "-"))

# Create the ggplot visualization
p2 <- ggplot(df_summary1, aes(x = Date_monthly, y = mean)) +
  geom_line(color = "#FC8D59", linewidth=2) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), fill = "#FC8D59", alpha = 0.5,color = "black", linetype = "dotted") +
  geom_line(aes(ymin = CI_lower, ymax = CI_upper)) + 
  geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),fill="white",alpha=1,inherit.aes = FALSE) +
  geom_vline(xintercept = rect$xmin, linetype="dashed", 
             color = "darkgrey") +
  geom_vline(xintercept = rect$xmax, linetype="dashed", 
             color = "darkgrey") +
  
  
  # Set x axis labels
  scale_x_date(limits = c(min(df_tidy$year_rev), max(df_tidy$year_rev)), date_breaks = "1 year", date_labels = "%Y")+ 
  
  ## Set the Y-axis scale, remove for autoscale
  
  coord_cartesian(ylim = c(-21, -15)) +
  
  scale_y_continuous(breaks=c(-21, -20, -19, -18, -17,-16,-15)) +
  ## Set theme&basic font
  theme_article(base_size = 20) +
  
  
  ### Style the axis (font size)
  # theme(axis.text.x = element_text(size=16, angle=0, vjust = 0), axis.text.y = element_text(size=16)) +
  
  ### Set layout of the graph 
  # theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black", fill=NA)) +
  
  ### Set label of x- and y-axis
  ylab(expression(paste(delta^{13}, "C (\u2030)"))) + xlab("Estimated date (years)") +
  
  ### Set aspect ratio of the graph n/n = square
  theme(aspect.ratio=2/4)

p2

# Create the ggplot visualization
#p2 <- ggplot(df_summary1, aes(x = Date_monthly, y = mean)) +
  geom_line() +
  #geom_line(data=df_tidy1, aes(x=as.Date(df_tidy1$year_rev), y=df_tidy1$d13c.cor, group=Whale), color="lightgrey", size=0.25, alpha=0.3) +
  geom_line(size = 0.9, alpha = 0.8, color="#FC8D59") +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), fill = "#FC8D59", color= "#FC8D59", linetype= "dotted", alpha = 0.3) +
  geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),fill="white",alpha=1,inherit.aes = FALSE) +
  geom_vline(xintercept = rect$xmin, linetype="dashed", 
             color = "darkgrey") +
  geom_vline(xintercept = rect$xmax, linetype="dashed", 
             color = "darkgrey") +
  
  
  # Set x axis labels
  scale_x_date(limits = c(min(df_tidy$year_rev), max(df_tidy$year_rev)), date_breaks = "1 year", date_labels = "%Y")+ 
  
  ## Set the Y-axis scale, remove for autoscale
  
  coord_cartesian(ylim = c(-21, -15)) +
  
  scale_y_continuous(breaks=c(-21, -20, -19, -18, -17,-16,-15)) +
  ## Set theme&basic font
  theme_article(base_size = 20) +
  
  
  ### Style the axis (font size)
  # theme(axis.text.x = element_text(size=16, angle=0, vjust = 0), axis.text.y = element_text(size=16)) +
  
  ### Set layout of the graph 
  # theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black", fill=NA)) +
  
  ### Set label of x- and y-axis
  ylab(expression(paste(delta^{13}, "C (\u2030)"))) + xlab("Year") +
  
  ### Set aspect ratio of the graph n/n = square
  theme(aspect.ratio=2/4)

p2

library(patchwork)
p1+ theme_article(base_size = 20)/p2+ theme_article(base_size = 20)

#--------------------------------------------------------------------------
#plot
library(ggplot2)
point_size=1
isop <- ggplot(fw.data_d13cor, aes(x = Cm)) + 
  theme_article(base_size = 20) + 
  xlab("Distance from gingiva (Cm)") + theme_article() + 
  scale_x_continuous(breaks=c(-10, 0, 10,20,30,40,50))


isop <- isop + geom_line(aes(y = dC),col = "#FC8D59") + 
  geom_point(aes(y = dC), size = point_size, shape = 16,col = "#FC8D59") +
  ylab(expression(paste(delta^{13}, "C (\u2030)"))) + facet_wrap(vars(Whale), ncol=8)

# add the d15N data but need to scale it to d13C
a <- mean(na.omit(fw.data_d13cor$dN))
a
b <- mean(na.omit(fw.data_d13cor$d13c.cor))
b

v_shift <- a - b

#v-shift (ind F13065) = v_shift

isop <- isop + geom_line(aes(y = dN - v_shift), col = "#74A9CF") + 
  geom_point(aes(y = (dN - v_shift)), size = point_size, 
             shape = 16, col = "#74A9CF")+ facet_wrap(vars(Whale), ncol = 6) + 
  scale_y_continuous(sec.axis = sec_axis(~.+v_shift, name = expression(paste(delta^{15}, "N (\u2030)")))) +
  
  ### Set aspect ratio of the graph n/n = square
  theme(aspect.ratio=2/4)
isop + theme_article(base_size = 20)
