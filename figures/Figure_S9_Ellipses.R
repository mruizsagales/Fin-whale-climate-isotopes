
#--------------------------------------------------------------------------------
# Figure S9 - Ellipses and SEAB/SEAC through time
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
library(coda)
library(tidyr)
library(tidybayes)
library(ellipse)

standard_theme <- theme(
  legend.text = element_text(size = 12),  
  legend.key.size = unit(0.5, 'cm'), 
  legend.key.height = unit(0.5, 'cm'), 
  legend.key.width = unit(0.5, 'cm'),
  legend.title = element_blank(),
  aspect.ratio=1,
  axis.text.x = element_text(angle = 45), 
  legend.position = "none",
  legend.direction = "vertical")


# 2. Import stable isotope data

df <- read_excel("~/Desktop/Doctorat/Analisis_isotops_barbes/Projecte_barbes_clima/Environmental_vars/All_merged_time_calibrated_2013_to_2022_Suess_cor.xlsx") # Import Suess-corrected dataset of stable isotope data along the baleen plate of fin whales 

# 3. Stable isotope ratios per whale and month

nb.cols <- 29
mycolors <- colorRampPalette(brewer.pal(8, "RdYlBu"))(nb.cols) #"YlGnBu"


a <- ggplot(df, aes(as.numeric(Month_from_sample_date), dN)) + theme_article(base_size = 15) + standard_theme + geom_smooth(aes(color=Whale), se = FALSE) + geom_smooth(color= "black", size=2.5) + scale_x_continuous(breaks = c(1:12), labels = c("Jan", "Feb","Mar","Apr", "May", "Jun","Jul","Aug","Sep","Oct","Nov","Dec")) + xlab("Month") + ylab(expression(paste(delta^{15}, "N (\u2030)"))) + scale_color_manual(values = mycolors)
b <- ggplot(df, aes(as.numeric(Month_from_sample_date), d13c.cor)) + theme_article(base_size = 15) + standard_theme + geom_smooth(aes(color=Whale), se = FALSE) + geom_smooth(color= "black", size=2.5)+ scale_x_continuous(breaks = c(1:12), labels = c("Jan", "Feb","Mar","Apr", "May", "Jun","Jul","Aug","Sep","Oct","Nov","Dec")) + xlab("Month") + ylab(expression(paste(delta^{13}, "C (\u2030)"))) + scale_color_manual(values = mycolors)

c<- a+b + patchwork::plot_annotation(tag_levels = "a", 
                                   tag_suffix = ")")
print(c)

ggsave("/Users/marcruizisagales/Documents/GitHub/Climate-baleen-plate-isotopes/png/Figure_7_Baleen_stable_isotope_ratios_for_each_whale_and_month.png", last_plot(), 
       dpi = 300,  width = 20, height = 15, units = "cm")

# 4. Isotopic range: overall and per whale
df1 <- df %>% dplyr::group_by(Whale) %>% dplyr::reframe(range_dN= max(dN, na.rm = TRUE)-min(dN, na.rm = TRUE), range_dC= max(d13c.cor, na.rm = TRUE)- min(d13c.cor, na.rm = TRUE))
df1

#overall
range(df$dN, na.rm = T)
range(df$d13c.cor, na.rm = T)

#per whale
mean(df1$range_dN, na.rm = T)
sd(df1$range_dN, na.rm = T)

mean(df1$range_dC, na.rm = T)
sd(df1$range_dC, na.rm = T)

# 4. Ellipses

standard_theme <- theme(
  legend.text = element_text(size = 12),  # Legend text
  legend.key.size = unit(0.5, 'cm'), 
  legend.key.height = unit(0.5, 'cm'), 
  legend.key.width = unit(0.5, 'cm'),
  legend.position = "bottom",
  legend.direction = "horizontal", 
  legend.title = element_blank())

# Clean dataset

df_1 <- df[complete.cases(df$dN),] # remove NA
df_2 <- df_1[complete.cases(df_1$d13c.cor),] # remove NA
df_2 <- subset(df_2, !(Year_from_sample_date == 2010 & Whale == "F13065") & !(Year_from_sample_date == 2010 & Whale == "F13076")& !(Year_from_sample_date == 2012 & Whale == "F15083")& !(Year_from_sample_date == 2012 & Whale == "F15097")) # remove due to low data of the whale that year
df_3 <- df_2[, c("d13c.cor", "dN", "Year_from_sample_date", "Species")]
names(df_3) <- c("iso1","iso2", "group", "community")
df_4 <- as.data.frame(df_3)

# create the siber object
siber.fin <- createSiberObject(df_4)

# Calculate summary statistics for each group: TA, SEA and SEAc
group.ML <- groupMetricsML(siber.fin)
print(group.ML)

df_group.ML <- as.data.frame(group.ML)
df_group.ML$RowNames <- rownames(df_group.ML)
reshaped_dataframe <- gather(df_group.ML, key = "Key", value = "Value", -RowNames) # Reshape the dataframe
reshaped_dataframe <- separate(reshaped_dataframe, Key, into = c("Whale", "Year"), sep = "\\.", remove = FALSE) # Separate the 'Key' column into 'Whale' and 'Year'

# Communuty metrics ML
community.ML <- communityMetricsML(siber.fin) 
print(community.ML)

# Plot ellipses

sbg_1 <- df %>% 
  dplyr::group_by(Year_from_sample_date, Whale) %>% 
  dplyr::summarise(mC = mean(d13c.cor, na.rm=TRUE), 
                   sdC = sd(d13c.cor, na.rm=TRUE), 
                   mN = mean(dN, na.rm=TRUE), 
                   sdN = sd(dN, na.rm=TRUE) )



coul <- brewer.pal(9, "RdYlBu") # set colors
coul <- colorRampPalette(coul)(30)
values <- setNames(coul, sort(unique(sbg_1$Whale)))

first.plot <- ggplot(data = as.data.frame(df), # plot points
                     mapping = aes(x = d13c.cor, 
                                   y = dN)) + 
  geom_point(aes(color = Whale), size = 0.22, alpha=0.9) +
  facet_wrap(~Year_from_sample_date, ncol = 4) +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) + theme_article(base_size = 15) +
  scale_colour_manual(values = values)

print(first.plot)

second.plot <- first.plot + # add errorbars
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

p.ell <- 0.95 # decide how big an ellipse you want to draw

ellipse.plot <- second.plot + # addellipses
  stat_ellipse(aes(group = Year_from_sample_date), alpha=0.01, 
               linetype = 2,
               color="darkgrey",
               level = p.ell,
               type = "norm",
               geom = "polygon")+ ylim(6,14)



ellipse.plot.final <- ellipse.plot + standard_theme + ### Set aspect ratio of the graph n/n = square
  theme(aspect.ratio=3/4) 

print(ellipse.plot.final) 

ggsave("/Users/marcruizisagales/Documents/GitHub/Climate-baleen-plate-isotopes/png/Figure_S9_Ellipses.png", last_plot(), 
       dpi = 300,  width = 20, height = 15, units = "cm")


# 5. SEAc in time

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

print(SEAc_time)

# 6. Mean dN vs SD dN

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
  xlab(expression(paste("Annual ",delta^{15}, "N (\u2030) mean"))) + theme_article(base_size = 15) + theme(aspect.ratio = 1)

print(Nmean_vs_NSD)

ggsave("/Users/marcruizisagales/Documents/GitHub/Climate-baleen-plate-isotopes/png/Figure_S9_Nmean_vs_NSD.png", last_plot(), 
       dpi = 300,  width = 20, height = 15, units = "cm")

# 7. SEAB

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

# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.

# get a list of all the files in the save directory
#all.files <- dir(parms$save.dir, full.names = TRUE)

# find which ones are jags model files
#model.files <- all.files[grep("jags_output", all.files)]

# test convergence for the first one
#do.this <- 1

#load(model.files[do.this])

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
write.xlsx(x=SEAB_posterior_estimations_for_each_year, file= "/Users/marcruizisagales/Documents/GitHub/Climate-baleen-plate-isotopes/data/SEAB_posterior_estimations_for_each_year.xlsx",overwrite = TRUE)

cr.p <- c(0.95, 0.99) # Calculate some credible intervals, vector of quantiles

SEA.B.credibles <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)

SEA.B.modes <- lapply(# do similar to get the modes, taking care to pick up multimodal posterior distributions if present
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = cr.p, all.modes=T)

mu.post <- extractPosteriorMeans(siber.fin, ellipses.posterior) # extract the posterior means

# calculate the corresponding distribution of layman metrics
layman.B <- bayesianLayman(mu.post)

#plot
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

df_long <- as.data.frame(df_ordered) %>%
  pivot_longer(cols = 1:13, names_to = "year", values_to = "value")

data_Seab <- df_long %>% dplyr::group_by(year) %>% dplyr::summarise(mean(value)) 


RColorBrewer::brewer.pal(5, "Reds")
SEAB_year <- ggplot() + stat_eye(data= df_long, aes(x = as.character(year), y = value), .width= c(0.5,0.75, 0.95), fill="#FCAE91")+ theme_article(base_size = 15) + geom_point(data= long_data1b, aes(x=as.character(long_data1b$year), y=value), color= "#FB6A4A") + theme(aspect.ratio = 3/4) + ylim(0,3.3) + xlab("Year") + ylab(expression("Standard Ellipse Area " ("\u2030"^2) )) + scale_x_discrete(labels= c(2010:2022))
print(SEAB_year)

ggsave("/Users/marcruizisagales/Documents/GitHub/Climate-baleen-plate-isotopes/png/Figure_7_Posterior_estimates_SEAB_per_year.png", last_plot(), 
       dpi = 300,  width = 20, height = 15, units = "cm")


# 8. GAM of SEAB per year

sbg_2 <- df_2 %>% dplyr::group_by(Year_from_sample_date) %>% dplyr::summarise(count = n(),
                                                                              mC = mean(d13c.cor), 
                                                                              sdC = sd(d13c.cor), 
                                                                              mN = mean(dN), 
                                                                              sdN = sd(dN))
colnames(sbg_2) <- c("Year","count","mC","sdC","mN","sdN")

data_Seab_df <- as.data.frame(data_Seab)
colnames(data_Seab_df) <- c("Year","mean(value)")

Seabb <- merge(sbg_2, data_Seab_df, by= "Year")

#as numeric
Seabb$Year <- as.double(Seabb$Year)
Seabb$Value <- as.numeric(Seabb$`mean(value)`)
Seabb$count <- as.numeric(Seabb$count)

#do gam
Seabb_gam <- mgcv::gam(Value ~  s(Year) + s(count), select=TRUE, method = 'REML',data= Seabb)
summary(Seabb_gam)

par(mar=c(5,5,5,5))
plot.gam(Seac_gam)

#predict
model_2_p <- predict_gam(Seabb_gam)

SEAB_time <- predict_gam(Seabb_gam, exclude_terms = "s(count)") %>%
  ggplot(aes(Year, fit)) + geom_smooth_ci() + geom_line(color = "darkgrey", linewidth=2) +
  geom_ribbon(aes(ymin = fit-1.96 *se.fit, ymax = fit+1.96 *se.fit), alpha = 0,color = "black", linetype = "dotted") + 
  geom_point(data = Seabb, aes(Year, Value)) +
  ylab(expression(paste("SEAB (\u2030"^2, ")" ))) +
  xlab("Year") + theme_article(base_size = 15) + theme(aspect.ratio = 1) + scale_x_continuous(breaks= c(2010, 2012, 2014,2016,2018,2020,2022))

print(SEAB_time)

ggsave("/Users/marcruizisagales/Documents/GitHub/Climate-baleen-plate-isotopes/png/Figure_S9_SEAB_time.png", last_plot(), 
       dpi = 300,  width = 20, height = 15, units = "cm")

# 9. Bayesian overlap among years

group_names <- names(ellipses.posterior) # Extract group names

overlap_matrix <- matrix(0, nrow = length(group_names), ncol = length(group_names), dimnames = list(group_names, group_names)) # Initialize a matrix to store overlap results (mean)

cred_matrix <- matrix(0, nrow = length(group_names), ncol = length(group_names), dimnames = list(group_names, group_names)) # Initialize a matrix to store overlap results (cred. int.)

#-
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

#-

colnames(overlap_df) <- gsub("fin.", "", colnames(overlap_df))
rownames(overlap_df) <- gsub("fin.", "", rownames(overlap_df))
order <- as.character(2010:2022)

overlap_df1 <- overlap_df[order, order]

print(overlap_df1)

cred_df <- as.data.frame(cred_matrix) 

colnames(cred_df) <- gsub("fin.", "", colnames(cred_df))
rownames(cred_df) <- gsub("fin.", "", rownames(cred_df))
order <- as.character(2010:2022)

cred_df1 <- cred_df[order, order]

print(cred_df1)


      
# 10. Draw ellipses

n.posts <- 10 # numer of posterior draws 
p.ell <- 0.95 # ellipsis size
#p.ell <- pchisq(1,2) # for a standard ellipse use

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
  geom_point(aes(color = factor(Year_from_sample_date):factor(Species)), size = 0.2)+
  ylab(expression(paste(delta^{15}, "N (permille)")))+
  xlab(expression(paste(delta^{13}, "C (permille)"))) + 
  theme(text = element_text(size=15), aspect.ratio = 1) + facet_wrap(vars(group))
print(first.plot)

second.plot <- first.plot + 
  geom_polygon(data = ellipse_df,
              mapping = aes(iso1, iso2,
                             group = group,
                             color = factor(group):factor(community),
                             fill = NULL),
               fill = NA,
               alpha = 1) + facet_wrap(vars(group)) + theme_article(base_size = 15) + theme(text = element_text(size=15), aspect.ratio = 1)
print(second.plot)


