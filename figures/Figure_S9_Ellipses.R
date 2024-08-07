
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


nb.cols <- 114
mycolors <- colorRampPalette(brewer.pal(8, "RdYlBu"))(nb.cols) #"YlGnBu"

df$Year_Whale <- paste0(df$Year_from_sample_date, sep="-", df$Whale)


a <- ggplot(df, aes(as.numeric(Month_from_sample_date), dN)) + theme_article(base_size = 15) + theme(aspect.ratio = 2/1,axis.text.x = element_text(angle = 45), legend.position = "none") + geom_smooth(aes(color=Year_Whale), se = FALSE) + geom_smooth(color= "black", size=2.5) + scale_x_continuous(breaks = c(1:12), labels = c("Jan", "Feb","Mar","Apr", "May", "Jun","Jul","Aug","Sep","Oct","Nov","Dec")) + xlab("Month") + ylab(expression(paste(delta^{15}, "N (\u2030)"))) + scale_color_manual(values = mycolors)
b <- ggplot(df, aes(as.numeric(Month_from_sample_date), d13c.cor)) + theme_article(base_size = 15) + theme(aspect.ratio = 2/1,axis.text.x = element_text(angle = 45), legend.position = "none") + geom_smooth(aes(color=Year_Whale), se = FALSE) + geom_smooth(color= "black", size=2.5)+ scale_x_continuous(breaks = c(1:12), labels = c("Jan", "Feb","Mar","Apr", "May", "Jun","Jul","Aug","Sep","Oct","Nov","Dec")) + xlab("Month") + ylab(expression(paste(delta^{13}, "C (\u2030)"))) + scale_color_manual(values = mycolors)

a/b
df1 <- df %>% dplyr::group_by(Whale) %>% dplyr::reframe(range_dN= max(dN, na.rm = TRUE)-min(dN, na.rm = TRUE), range_dC= max(d13c.cor, na.rm = TRUE)- min(d13c.cor, na.rm = TRUE))

ggplot(df, aes(Month_from_sample_date, dN)) + geom_boxplot()

df %>% dplyr::group_by(Month_from_sample_date) %>% dplyr::summarise(mean_dN=mean(dN, na.rm=T), mean_d13c.cor=mean(d13c.cor, na.rm=T))
10.3-9.17
-19.3-(-19.6)

range(df$dN, na.rm = T)
range(df$d13c.cor, na.rm = T)

mean(df1$range_dN, na.rm = T)
sd(df1$range_dN, na.rm = T)

mean(df1$range_dC, na.rm = T)
sd(df1$range_dC, na.rm = T)

mean(c(1,0.5,1.6,0.9,1,1.4,1.3,0.9,1.5,2.2))
sd(c(1,0.5,1.6,0.9,1,1.4,1.3,0.9,1.5,2.2))
length(c(1,0.5,1.6,0.9,1,1.4,1.3,0.9,1.5,2.2))

mean(c(0.8,1,1,1.1,0.8,1,0.6,1.3,0.8,1.9))
sd(c(0.8,1,1,1.1,0.8,1,0.6,1.3,0.8,1.9))
length(c(0.8,1,1,1.1,0.8,1,0.6,1.3,0.8,1.9))

DataDryad <- read.csv("~/Downloads/DataDryad.csv", header=FALSE)

DataDryad %>% filter(Species == "Balaenoptera physalus") %>% dplyr::group_by(Lab.ID) %>% dplyr::reframe(rN=range(max(ẟ15N)-min(ẟ15N), na.rm=TRUE), rC=range(max(ẟ13C)-min(ẟ13C)), na.rm=TRUE)

df_1 <- df[complete.cases(df$dN),] # remove NA
df_2 <- df_1[complete.cases(df_1$d13c.cor),] # remove NA

sort(unique(df_2$Year_from_sample_date)) # years

mean(c(3.11,2.29,3.57,2.95,2.92))
sd(c(3.11,2.29,3.57,2.95,2.92))
length(c(3.11,2.29,3.57,2.95,2.92))

mean(c(5.7,6.83,6.33,5.11,8.59))
sd(c(5.7,6.83,6.33,5.11,8.59))
length(c(5.7,6.83,6.33,5.11,8.59))


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



ellipse.plot.final <- ellipse.plot + theme(legend.position = c(0.6, 0.1),
                            legend.direction = "horizontal", legend.title = element_blank()) + ### Set aspect ratio of the graph n/n = square
  theme(aspect.ratio=2/4) 

print(ellipse.plot.final) 

ggsave("/Users/marcruizisagales/Documents/GitHub/Climate-baleen-plate-isotopes/png/Figure_S9_Ellipses.png", ellipse.plot.final, 
       device = png(width = 1200, height = 450))


# 2. Plot 2

sbg_2 <- df_2 %>% 
  dplyr::group_by(Year_from_sample_date) %>% 
  dplyr::summarise(count = n(),
                   mC = mean(d13c.cor), 
                   sdC = sd(d13c.cor), 
                   mN = mean(dN), 
                   sdN = sd(dN))

colnames(sbg_2) <- c("Year","count","mC","sdC","mN","sdN")
df_SEAc <- filter(reshaped_dataframe, RowNames == "SEAc")

Seac <- merge(sbg_2, df_SEAc, by= "Year")



Seac$Year <- as.double(Seac$Year)
Seac$Value <- as.numeric(Seac$Value)
Seac$count <- as.numeric(Seac$count)

Seac_gam <- mgcv::gam(Value ~  s(Year) + s(count), select=TRUE, method = 'REML',
                        data= Seac)

summary(Seac_gam)

par(mar=c(5,5,5,5))
plot.gam(Seac_gam)

model_2_p <- predict_gam(Seac_gam)
model_2_p


SEAc_time <- predict_gam(Seac_gam, exclude_terms = "s(count)") %>%
  ggplot(aes(Year, fit)) + geom_smooth_ci() + geom_line(color = "darkgrey", linewidth=2) +
  geom_ribbon(aes(ymin = fit-1.96 *se.fit, ymax = fit+1.96 *se.fit), alpha = 0,color = "black", linetype = "dotted") + 
  geom_point(data = Seac, aes(Year, Value)) +
  ylab(expression(paste("SEAc (\u2030"^2, ")" ))) +
  xlab("Estimated date (years)") + theme_article(base_size = 15) + theme(aspect.ratio = 1) + scale_x_continuous(breaks= c(2010, 2012, 2014,2016,2018,2020,2022))

ggsave("/Users/marcruizisagales/Documents/GitHub/Climate-baleen-plate-isotopes/png/Figure_S9_SEAc_time.png", SEAc_time, 
       device = png(width = 1200, height = 450))

# 3. Plot 3

sbg_2$Year <- as.numeric(sbg_2$Year)
sbg_2$sdN <- as.numeric(sbg_2$sdN)
sbg_2$mN <- as.numeric(sbg_2$mN)

SDn_mN_gam <- mgcv::gam(sdN ~  s(mN) + s(count),
                                data= sbg_2,method="REML", select=TRUE)

summary(SDn_mN_gam)
par(mar=c(5,5,5,5))
plot.gam(SDn_mN_gam)

model_2_p <- predict_gam(SDn_mN_gam)
model_2_p


Nmean_vs_NSD <- predict_gam(SDn_mN_gam, exclude_terms = "s(count)") %>%
  ggplot(aes(mN, fit)) + geom_smooth_ci() + geom_line(color = "darkgrey", linewidth=2) +
  geom_ribbon(aes(ymin = fit-1.96 *se.fit, ymax = fit+1.96 *se.fit), alpha = 0,color = "black", linetype = "dotted") + 
  ylab(expression(paste("Annual ",delta^{15}, "N (\u2030) SD"))) +
  xlab(expression(paste("Annual ",delta^{15}, "N (\u2030) mean"))) + theme_article(base_size = 20) + theme_article(base_size = 15) + theme(aspect.ratio = 1)

ggsave("/Users/marcruizisagales/Documents/GitHub/Climate-baleen-plate-isotopes/png/Figure_S9_Nmean_vs_NSD.png", Nmean_vs_NSD, 
       device = png(width = 1200, height = 450))

#-------
# SEAb
#-------

# options for running jags
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.
ellipses.posterior <- siberMVN(siber.fin, parms, priors)
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 30
#>    Unobserved stochastic nodes: 3
#>    Total graph size: 45
#> 
#> Initializing model
#> 
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 30
#>    Unobserved stochastic nodes: 3
#>    Total graph size: 45
#> 
#> Initializing model
#> 
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 30
#>    Unobserved stochastic nodes: 3
#>    Total graph size: 45
#> 
#> Initializing model
#> 
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 30
#>    Unobserved stochastic nodes: 3
#>    Total graph size: 45
#> 
#> Initializing model
#> 
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 30
#>    Unobserved stochastic nodes: 3
#>    Total graph size: 45
#> 
#> Initializing model
#> 
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 30
#>    Unobserved stochastic nodes: 3
#>    Total graph size: 45
#> 
#> Initializing model
#> 
#> # The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.

# get a list of all the files in the save directory
all.files <- dir(parms$save.dir, full.names = TRUE)

# find which ones are jags model files
model.files <- all.files[grep("jags_output", all.files)]

# test convergence for the first one
do.this <- 1

load(model.files[do.this])

gelman.diag(output, multivariate = FALSE)

gelman.plot(output, auto.layout = FALSE)

SEA.B <- siberEllipses(ellipses.posterior)
colnames(SEA.B) <- colnames(group.ML) 

years <- as.numeric(sub("fin\\.", "", colnames(SEA.B))) # Extract the years from the column names

# Order the columns based on the extracted years
ordered_cols <- colnames(SEA.B)[order(years)]

# Reorder the dataframe columns
df_ordered <- SEA.B[, ordered_cols]
df_ordered1 <- group.ML[, ordered_cols]
colnames(df_ordered) <- gsub("fin.", "", colnames(df_ordered))
SEAB_posterior_estimations_for_each_year <- as.data.frame(df_ordered)

library(openxlsx)
write.xlsx(x=SEAB_posterior_estimations_for_each_year, file= "/Users/marcruizisagales/Desktop/SEAB_posterior_estimations_for_each_year.xlsx",overwrite = TRUE)

siberDensityPlot(df_ordered, xticklabels = colnames(df_ordered), 
                 xlab = c("Year"),
                 ylab = expression("Standard Ellipse Area " ("\u2030"^2) ),
                 bty = "o",
                 las = 1,
                 ylims =c(0,3.3)
)

# Add red x's for the ML estimated SEA-c
points(1:ncol(df_ordered), df_ordered1[3,], col="red", pch = "x", lwd = 2)

df_ordered1
colnames(df_ordered1) <- gsub("fin.", "", colnames(df_ordered1))

# Create the initial dataframe
data <- data.frame(
  year = c(2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022),
  TA = c(6.366250, 9.232935, 6.300406, 7.911607, 8.532122, 8.559046, 15.997114, 7.027397, 5.348467, 1.5765905, 3.9680961, 3.7752765, 4.5543600),
  SEA = c(1.936771, 2.236560, 1.006723, 1.330005, 1.313588, 1.638420, 2.615495, 1.384236, 1.144628, 0.5458290, 0.9056720, 0.9701721, 0.9053202),
  SEAc = c(2.005941, 2.260871, 1.013813, 1.336897, 1.322904, 1.647946, 2.631944, 1.393052, 1.158587, 0.5628862, 0.9205191, 0.9858200, 0.9319473)
)

# Reshape the dataframe from wide to long format
long_data1 <- pivot_longer(data, cols = -year, names_to = "class", values_to = "value")

# Display the result
print(long_data1)

long_data1b <- long_data1 %>% filter(class == "SEAc")
#
library(tidyr)
df_long <- as.data.frame(df_ordered) %>%
  pivot_longer(cols = 1:13, names_to = "year", values_to = "value")
library(tidybayes)
RColorBrewer::brewer.pal(5, "Reds")
ggplot() +
  stat_eye(data= df_long, aes(x = as.character(year), y = value), .width= c(0.5,0.75, 0.95), fill="#FCAE91")+ theme_article(base_size = 15) + geom_point(data= long_data1b, aes(x=as.character(long_data1b$year), y=value), color= "#FB6A4A") + theme(aspect.ratio = 3/4) + ylim(0,3.3) + xlab("Year") + ylab(expression("Standard Ellipse Area " ("\u2030"^2) )) + scale_x_discrete(labels= c(2010:2022))

data_Seab <- df_long %>% dplyr::group_by(year) %>% dplyr::summarise(mean(value)) 

# Calculate some credible intervals 
cr.p <- c(0.95, 0.99) # vector of quantiles

# call to hdrcde:hdr using lapply()
SEA.B.credibles <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)

# do similar to get the modes, taking care to pick up multimodal posterior
# distributions if present
SEA.B.modes <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = cr.p, all.modes=T)

# extract the posterior means
mu.post <- extractPosteriorMeans(siber.example, ellipses.posterior)

# calculate the corresponding distribution of layman metrics
layman.B <- bayesianLayman(mu.post)

### Bayesian
#example
bayes.overlap.G2.G3 <- bayesianOverlap("fin.2011", "fin.2012", ellipses.posterior, 
draws = 10, p.interval = 0.95,
n = 360)
print(bayes.overlap.G2.G3)

mean(bayes.overlap.G2.G3$overlap)

bayes.prop.95.over <- (bayes.overlap.G2.G3[,3] / (bayes.overlap.G2.G3[,2] + 
                                                    bayes.overlap.G2.G3[,1] -
                                                    bayes.overlap.G2.G3[,3])*100
)

hist(bayes.prop.95.over, 10)
mean(bayes.prop.95.over)

SEA.B.credibles <- lapply(
  as.data.frame(bayes.prop.95.over), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)

library(tidyr)




# Extract group names
group_names <- names(ellipses.posterior)

# Initialize a matrix to store overlap results
overlap_matrix <- matrix(0, nrow = length(group_names), ncol = length(group_names), 
                         dimnames = list(group_names, group_names))

cred_matrix <- matrix(0, nrow = length(group_names), ncol = length(group_names), 
                         dimnames = list(group_names, group_names))


# Loop to calculate overlaps for all pairs
for (i in 1:(length(group_names))) {
  for (j in 1:(length(group_names))) {
    a <- bayesianOverlap(group_names[i], group_names[j], ellipses.posterior, draws = 10, p.interval = 0.95,
                         n = 100)
    
    bayes.prop.95.over <- (a[,3] / (a[,2] + 
                                      a[,1] -
                                      a[,3])*100
    )
    
    mean_a <- mean(bayes.prop.95.over)
    
    SEA.B.credibles <- lapply(
      as.data.frame(bayes.prop.95.over), 
      function(x,...){tmp<-hdrcde::hdr(x)$hdr},
      prob = cr.p)
    
    cred_a<-  paste(round(SEA.B.credibles$bayes.prop.95.over[2,1],2), round(SEA.B.credibles$bayes.prop.95.over[2,2],2), sep = "-")
    
    overlap_matrix[i, j] <- mean_a
    overlap_matrix[j, i] <- overlap_matrix[i, j]  # Symmetric matrix
    
    cred_matrix[i, j] <- cred_a
    cred_matrix[j, i] <- cred_matrix[i, j]  # Symmetric matrix
  }
}

# Convert matrix to data frame for better readability
overlap_df <- as.data.frame(overlap_matrix)

# Display the overlap table
print(overlap_df)

# Remover el prefijo "fin." de los nombres de columnas
colnames(overlap_df) <- gsub("fin.", "", colnames(overlap_df))

# Remover el prefijo "fin." de los nombres de filas
rownames(overlap_df) <- gsub("fin.", "", rownames(overlap_df))

print(overlap_df)

# Especificar el orden deseado
order <- as.character(2010:2022)

# Reordenar las columnas y filas
overlap_df1 <- overlap_df[order, order]

print(overlap_df1)


# Convert matrix to data frame for better readability
cred_df <- as.data.frame(cred_matrix)

# Display the overlap table
print(cred_df)

# Remover el prefijo "fin." de los nombres de columnas
colnames(cred_df) <- gsub("fin.", "", colnames(cred_df))

# Remover el prefijo "fin." de los nombres de filas
rownames(cred_df) <- gsub("fin.", "", rownames(cred_df))

print(cred_df)

# Especificar el orden deseado
order <- as.character(2010:2022)

# Reordenar las columnas y filas
cred_df1 <- cred_df[order, order]

print(cred_df1)


# a histogram of the overlap
hist(bayes.overlap.G2.G3[,3], 10)

# and as above, you can express this a proportion of the non-overlapping area of 
# the two ellipses, would be
bayes.prop.95.over <- (bayes.overlap.G2.G3[,3] / (bayes.overlap.G2.G3[,2] + 
                                                    bayes.overlap.G2.G3[,1] -
                                                    bayes.overlap.G2.G3[,3])
)

hist(bayes.prop.95.over, 10)
#-----------------------
      
# drow ellipses
library(ellipse)
# how many of the posterior draws do you want?
n.posts <- 10

# decide how big an ellipse you want to draw
p.ell <- 0.95

# for a standard ellipse use
# p.ell <- pchisq(1,2)




# a list to store the results
all_ellipses <- list()

# loop over groups
for (i in 1:length(ellipses.posterior)){

  # a dummy variable to build in the loop
  ell <- NULL
  post.id <- NULL

  for ( j in 1:n.posts){

    # covariance matrix
    Sigma  <- matrix(ellipses.posterior[[i]][j,1:4], 2, 2)

    # mean
    mu     <- ellipses.posterior[[i]][j,5:6]

    # ellipse points

    out <- ellipse::ellipse(Sigma, centre = mu , level = p.ell)


    ell <- rbind(ell, out)
    post.id <- c(post.id, rep(j, nrow(out)))

  }
  ell <- as.data.frame(ell)
  ell$rep <- post.id
  all_ellipses[[i]] <- ell
}

ellipse_df <- bind_rows(all_ellipses, .id = "id")


# now we need the group and community names

# extract them from the ellipses.posterior list
group_comm_names <- names(ellipses.posterior)[as.numeric(ellipse_df$id)]

# split them and conver to a matrix, NB byrow = T
split_group_comm <- matrix(unlist(strsplit(group_comm_names, "[.]")),
                           nrow(ellipse_df), 2, byrow = TRUE)

ellipse_df$community <- split_group_comm[,1]
ellipse_df$group     <- split_group_comm[,2]

ellipse_df <- dplyr::rename(ellipse_df, iso1 = x, iso2 = y)

df$group <- df$Year_from_sample_date
first.plot <- ggplot(data = as.data.frame(df), aes(d13c.cor, dN)) +
  geom_point(aes(color = factor(Year_from_sample_date):factor(Species)), size = 2)+
  ylab(expression(paste(delta^{15}, "N (permille)")))+
  xlab(expression(paste(delta^{13}, "C (permille)"))) + 
  theme(text = element_text(size=15)) + facet_wrap(vars(group))
print(first.plot)

second.plot <- first.plot + facet_wrap(~factor(Year_from_sample_date):factor(Species))
print(second.plot)

# rename columns of ellipse_df to match the aesthetics

third.plot <- first.plot + 
  geom_polygon(data = ellipse_df,
              mapping = aes(iso1, iso2,
                             group = group,
                             color = factor(group):factor(community),
                             fill = NULL),
               fill = NA,
               alpha = 1) + facet_wrap(vars(group)) + theme_article()
print(third.plot)

# 2. Plot 2

sbg_2 <- df_2 %>% 
  dplyr::group_by(Year_from_sample_date) %>% 
  dplyr::summarise(count = n(),
                   mC = mean(d13c.cor), 
                   sdC = sd(d13c.cor), 
                   mN = mean(dN), 
                   sdN = sd(dN))

colnames(sbg_2) <- c("Year","count","mC","sdC","mN","sdN")
data_Seab_df <- as.data.frame(data_Seab)
colnames(data_Seab_df) <- c("Year","mean(value)")

data_Seab_df$Year <- gsub("fin.", "", data_Seab_df$Year)

Seabb <- merge(sbg_2, data_Seab_df, by= "Year")



Seabb$Year <- as.double(Seabb$Year)
Seabb$Value <- as.numeric(Seabb$`mean(value)`)
Seabb$count <- as.numeric(Seabb$count)

Seabb_gam <- mgcv::gam(Value ~  s(Year) + s(count), select=TRUE, method = 'REML',
                      data= Seabb)

summary(Seabb_gam)

par(mar=c(5,5,5,5))
plot.gam(Seac_gam)

model_2_p <- predict_gam(Seabb_gam)
model_2_p


SEAc_time <- predict_gam(Seabb_gam, exclude_terms = "s(count)") %>%
  ggplot(aes(Year, fit)) + geom_smooth_ci() + geom_line(color = "darkgrey", linewidth=2) +
  geom_ribbon(aes(ymin = fit-1.96 *se.fit, ymax = fit+1.96 *se.fit), alpha = 0,color = "black", linetype = "dotted") + 
  geom_point(data = Seabb, aes(Year, Value)) +
  ylab(expression(paste("SEAB (\u2030"^2, ")" ))) +
  xlab("Year") + theme_article(base_size = 15) + theme(aspect.ratio = 1) + scale_x_continuous(breaks= c(2010, 2012, 2014,2016,2018,2020,2022))

