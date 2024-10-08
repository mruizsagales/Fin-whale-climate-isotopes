#--------------------------------------------------------------------------------
#Cross-correlation functions (CCF) in stable isotope oscillations in baleen plates
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
library(stats)
library(corrplot)
library(viridis)
library(ggcorrplot)

# 2. Import data

df <- read_excel("/Users/marcruizisagales/Documents/GitHub/Climate-baleen-plate-isotopes/data/0_point_alignments_29_juny_total_-4.xlsx") # Import dataset of the one-centimetre-spaced stable isotope data along the baleen plate of fin whales

min_cm_minus_4 <- df %>% dplyr::group_by(Whale) %>% dplyr::summarise(min_cm = min(Cm)) # Make sure that all the stable isotope data start at the position -4 cm 

# 3. Correlation between nitrogen stable isotope values in baleen from different whales

whale_names <- unique(df$Whale) # Get unique whale names

num_whales <- length(whale_names) # Number of individuals

correlation_matrix <- matrix(1,nrow = num_whales, ncol = num_whales, dimnames = list(whale_names, whale_names)) # Create an empty matrix to store max correlation values

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

# Plot the correlation matrix
corrplot(correlation_matrix, method="color", tl.col="black")

corplo_N <- ggcorrplot(correlation_matrix)
print(corplo_N)

ggsave("/Users/marcruizisagales/Documents/GitHub/Climate-baleen-plate-isotopes/png/Figure_S1_Corrplot_nitrogen.png", last_plot(), 
       dpi = 300,  width = 20, height = 20, units = "cm")

mean(correlation_matrix) # Average correlation between the nitrogen oscillations of different individuals

# 4. Correlation between carbon stable isotope values in baleen from different whales

whale_names <- unique(df$Whale) # Get unique whale names

num_whales <- length(whale_names) # Number of individuals

correlation_matrix <- matrix(1,nrow = num_whales, ncol = num_whales, dimnames = list(whale_names, whale_names)) # Create an empty matrix to store max correlation values

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

# Plot the correlation matrix
corrplot(correlation_matrix, method="color", tl.col="black")

corplo_C <- ggcorrplot(correlation_matrix)
print(corplo_C)

ggsave("/Users/marcruizisagales/Documents/GitHub/Climate-baleen-plate-isotopes/png/Figure_S2_Corrplot_carbon.png", last_plot(), 
       dpi = 300,  width = 20, height = 20, units = "cm")

mean(correlation_matrix) # Average correlation between the carbon oscillations of different individuals

# 5. Correlation between the nitrogen and carbon oscillations of the same individual

whale_names <- unique(df$Whale) # Get unique whale names

num_whales <- length(whale_names) # Number of individuals

correlation_matrix <- matrix(NA,nrow = num_whales, dimnames = list(whale_names)) # Create an empty matrix to store max correlation values

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

mean(correlation_matrix) # Average correlation between the nitrogen and carbon oscillations of the same individual

whales <- rownames(correlation_matrix) # whale code
ccf <- as.vector(correlation_matrix) # correlation N vs C
ccf_whales <- data.frame(whales, ccf) # dataframe of the N vs C correlation for each individual
range(ccf_whales$ccf) # range of the correlations

# plot of the N vs C correlation for each individual

NvC <- ggplot(ccf_whales, aes(x = whales, y = ccf, fill = ccf)) +
  geom_bar(stat = "identity") +
  scale_fill_distiller(palette = "RdBu", direction = 1) +
  scale_y_continuous(limits=c(0,1)) +
  labs(x = "Whales", y="") +
  theme_article(base_size = 15) +
  theme(axis.text.x=element_text(angle=90, hjust=1), legend.position = c(0.9,0.8), aspect.ratio = 3/4) + labs(fill = "CCF")

print(NvC)

ggsave("/Users/marcruizisagales/Documents/GitHub/Climate-baleen-plate-isotopes/png/Figure_S3_corr_N_C.png", NvC, 
       dpi = 300,  width = 20, height = 20, units = "cm")
