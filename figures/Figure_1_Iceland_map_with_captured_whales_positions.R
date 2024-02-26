
#--------------------------------------------------------------------------------
#Figure 1 - Iceland map with fin whale catch positions
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
library(vmsbase)
library(oce)
library(ocedata)
library(ncdf4)
library(RNetCDF)
library(sf)
library(rnaturalearth)
library(RColorBrewer)
library(marmap)

# 2. Import data

data <- read_excel("/Users/marcruizisagales/Desktop/Doctorat/Analisis_isotops_barbes/Projecte_barbes_clima/Inventari_barbes_Marc-2.xlsx")

data<- dplyr::filter(data, data$Feta == "Si") # sampled yes
data <- dplyr::filter(data, data$Hyb == "No") # hybrids no


# 3. Convert longitude and latitude from degrees, minutes and seconds to decimal units with ´lonsex2dec´ and ´latsex2dec´functions.

latsex2dec <- function(degree, minute, second, direction) {
  if (direction == "N") {
    declat <- degree + (minute / 60) + (second / 3600)
  } else {
    declat <- -(degree + (minute / 60) + (second / 3600))
  }
  return(declat)
}

# decimal latitude
data$LatDEC <- latsex2dec(degree = as.numeric(substr(data$Lat, 1, 2)), minute = as.numeric(substr(data$Lat, 3, 4)), second = 0, direction = "N")

lonsex2dec <- function(degree, minute, second, direction) {
  if (direction == "E") {
    declon <- degree + (minute / 60) + (second / 3600)
  } else {
    declon <- -(degree + (minute / 60) + (second / 3600))
  }
  return(declon)
}

#decimal longitude
data$LonDEC <- lonsex2dec(degree = as.numeric(substr(data$Long, 1, 2)), minute = as.numeric(substr(data$Long, 3, 4)), second = 0, direction = "W")


# 4. Create plot

points_labels <- # geographical labels
  data.frame( 
  lon = c(64.7,69, 62), 
  lat = c(-18,-35,-35.5),
  label = c("Iceland", "Greenland","North Atlantic"),
  face = c("plain", "plain", "italic"),
  color= c("black", "black", "black"),
  stringsAsFactors = FALSE
)

color= c("black", "black", "black") #for geographical labels

points_whales <- data.frame(lon = data$LonDEC, lat = data$LatDEC, var = c("2013", "2013","2013","2013","2013","2013","2015","2015","2015","2015","2015","2015","2015","2015","2015","2018","2018","2018","2018","2018","2018","2018","2018","2018","2018","2022","2022","2022","2022"))
points_whales$Whales <- data$ID

world <- ne_countries(scale = "large", returnclass = "sf") # world map

bathy_fname <- # Get bathymetric data
  "/Users/marcruizisagales/Downloads/GEBCO_20_Feb_2024_1040b342cfcb/gebco_2023_n70.0_s60.0_w-40.0_e-10.0.nc"
nc <- open.nc(bathy_fname)
tmp <- read.nc(nc)
tmp<- marmap::readGEBCO.bathy(bathy_fname)
#tmp <- getNOAA.bathy(-40, -15, 60, 70, res = 4, keep = TRUE)
bat_xyz <- as.xyz(tmp)
bat_xyz <- dplyr::filter(bat_xyz, V3<=0)

nb.cols <- 29
mycolors <- colorRampPalette(brewer.pal(9, "RdYlBu"))(nb.cols) # color palette points

nb.cols <- 100
a <- colorRampPalette(brewer.pal(9, "Blues"))(nb.cols) # color palette bathymetry

#alternativa_1
bat_labels <- data.frame(
  lon = c(-30), 
  lat = c(61),
  label = c("1000 m"),  # Replace with your labels
  stringsAsFactors = FALSE
)

b <- ggplot(data = world) +
  geom_tile(data = bat_xyz, aes(x = V1, y = V2, fill = V3), alpha=0.8) +
  #geom_tile(data = bat_xyz, aes(x = V1, y = V2, fill = V3), alpha=0.75) +
  geom_sf(data = world,fill = "grey", color= "black") +
  scale_fill_gradient("Depth (m)", low = a[99],high = a[1])+ 
  #geom_contour(data = bat_xyz, aes(x = V1, y = V2, z = V3),binwidth = 100, color = "grey85", size = 0.1) +
  #geom_contour(data = bat_xyz, aes(x = V1, y = V2, z = V3),breaks = -200, color = "grey85", size = 0.5) +
  coord_sf(xlim = c(-40, -10), ylim = c(60, 70),expand = FALSE) +
  #geom_point(data = points_labels, aes(x = lat, y = lon), color = "red", size = 1) +
  geom_text(data = points_labels, aes(x = lat, y = lon, label = label, fontface=face), vjust = -0.5, color=color) + labs(x = "Longitude", y = "Latitude", fill = "Whale ID") +
  geom_point(data=points_whales, aes(lon, lat), shape= 21, size=2, fill= mycolors, color="black")  + xlab("Longitude") + ylab("Latitude") +
  theme_article(base_size = 15) + theme(#aspect.ratio = 1,
                                        legend.position = c(0.9,0.1),
                                        legend.background = element_rect(fill=alpha('white', 0.2)),
                                        legend.text=element_text(size=8),
                                        legend.title=element_text(size=8),
                                        legend.key.size = unit(0.1, 'cm'), 
                                        legend.key.height = unit(0.2, 'cm'), 
                                        legend.key.width = unit(0.3, 'cm'),
                                        plot.margin = margin(2, 2, 2, 2, "cm"))#+ geom_text(data=bat_labels, aes(x = lon, y = lat, label = label), size=2.5, col="black")

b

ggsave("/Users/marcruizisagales/Documents/GitHub/Climate-baleen-plate-isotopes/png/Figure_1_Iceland_map_with_fin_whale_catch_positions.png", b, 
       device = png(width = 600, height = 400)
)


b <- ggplot(data = world) +
  #geom_tile(data = bat_xyz, aes(x = V1, y = V2, fill = V3), alpha=0.8) +
  #geom_tile(data = bat_xyz, aes(x = V1, y = V2, fill = V3), alpha=0.75) +
  geom_sf(data = world,fill = "grey", color= "black") +
  #scale_fill_gradient("Depth (m)", low = a[99],high = a[1])+ 
  geom_contour(data = bat_xyz, aes(x = V1, y = V2, z = V3),binwidth = 100, color = "grey85", linewidth = 0.1) +
  geom_contour(data = bat_xyz, aes(x = V1, y = V2, z = V3),breaks = -200, color = "grey85", linewidth = 0.5) +
  coord_sf(xlim = c(-40, -10), ylim = c(60, 70),expand = FALSE) +
  #geom_point(data = points_labels, aes(x = lat, y = lon), color = "red", size = 1) +
  geom_text(data = points_labels, aes(x = lat, y = lon, label = label, fontface=face), vjust = -0.5, color=color) + labs(x = "Longitude", y = "Latitude", fill = "Whale ID") +
  geom_point(data=points_whales, aes(lon, lat, fill=Whales), shape= 21, size=2)  + xlab("Longitude") + ylab("Latitude") +
  theme_article(base_size = 15) + theme(legend.background = element_rect(fill=alpha('white', 0.2)),
                                        legend.text=element_text(size=8),
                                        legend.title=element_text(size=8),
                                        legend.key.size = unit(0.1, 'cm'), #change legend key size
                                        legend.key.height = unit(0.2, 'cm'), #change legend key height
                                        legend.key.width = unit(0.3, 'cm'),
                                        plot.margin = margin(2, 2, 2, 2, "cm")) + scale_fill_manual("Whale ID", values=mycolors) + theme_article(base_size = 15) #+ geom_text(data=bat_labels, aes(x = lon, y = lat, label = label), size=2.5, col="black")

b

ggsave("/Users/marcruizisagales/Documents/GitHub/Climate-baleen-plate-isotopes/png/Figure_1_Iceland_map_with_fin_whale_catch_positions.png", b, 
       device = png(width = 600, height = 400)
)



