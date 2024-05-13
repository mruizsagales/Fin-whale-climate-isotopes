
#--------------------------------------------------------------------------------
# Figure 5 - Stable isotope of Northern krill from North Atlantic fin whales stomach contents
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
library(nortest)
library(car)
library(SuessR)

# 2. Import data

data <- read_excel("~/Desktop/Doctorat/Analisis_isotops_barbes/Clima/Data/Krill_Islandia.xlsx", sheet = "Hoja2")

mean(data$dN, na.rm=TRUE)
sd(data$dN, na.rm=TRUE)
range(data$dN, na.rm=TRUE)

mean(data$dC, na.rm=TRUE)
sd(data$dC, na.rm=TRUE)
range(data$dC, na.rm=TRUE)

# 3. Boxplot for nitrogen

krill_N <- ggplot(data, aes(as.factor(any),dN)) + 
  geom_boxplot(fill="#74A9CF", alpha=0.5, width=0.6) + 
  ylim(4,9) +
  geom_jitter(color="black", size=1, alpha=0.9, shape=1, width = 0.1) +
  #stat_summary(fun.y=mean, geom="point", shape=20, size=2, color="red", fill="red") +
  xlab("Year") + 
  ylab(expression(paste(delta^{15}, "N (‰)"))) + 
  theme_article(base_size = 15) + theme(aspect.ratio = 3/4)
krill_N

# 4. Boxplot for carbon
krill_C<- ggplot(data, aes(as.factor(any),dC)) + 
  geom_boxplot(fill="#74A9CF", alpha=0.5, width=0.6) +
  geom_jitter(color="black", size=1, alpha=0.9, shape=1, width = 0.1) +
  #stat_summary(fun.y=mean, geom="point", shape=20, size=2, color="red", fill="red") +
  xlab("Year") + 
  ylab(expression(paste(delta^{15}, "C (‰)"))) + 
  theme_article(base_size = 20) 
krill_C

#4.2. Boxplot for corrected carbon
data_cor <- data.frame(id = seq(1, length(data$dC)), year = data$any, region = "Subpolar North Atlantic", d13c = data$dC)
data_cor_suess <- SuessR(data_cor, correct.to = 2022)

mean(data_cor_suess$d13c.cor, na.rm=TRUE)
sd(data_cor_suess$d13c.cor, na.rm=TRUE)
range(data_cor_suess$d13c.cor, na.rm=TRUE)

krill_C_cor<- ggplot(data_cor_suess, aes(as.factor(year),d13c.cor)) + 
  geom_boxplot(fill="#74A9CF", alpha=0.5, width=0.6) +
  geom_jitter(color="black", size=1, alpha=0.9, shape=1, width = 0.1) +
  ylim(-21.5,-19) +
  #stat_summary(fun.y=mean, geom="point", shape=20, size=2, color="red", fill="red") +
  xlab("Year") + 
  ylab(expression(paste(delta^{15}, "C (‰)"))) + 
  theme_article(base_size = 15) 
krill_C_cor

plot <- krill_N + krill_C_cor +
  patchwork::plot_annotation(tag_levels = "a", 
                  tag_suffix = ")")

plot
ggsave("/Users/marcruizisagales/Documents/GitHub/Climate-baleen-plate-isotopes/png/Figure_5_krill_boxplots.png", plot, 
       device = png(width = 600, height = 400)
)

# 5. Perform ANOVA
data$any <- as.factor(data$any)
anova_result <- aov(dN ~ any, data = data)
summary(anova_result)

data_cor_suess$year <- as.factor(data_cor_suess$year)
anova_result <- aov(d13c.cor ~ year, data = data_cor_suess)
summary(anova_result)

# 6. Validation (Check assumptions)
# 6.1. Normality of Residuals
qqnorm(resid(anova_result))
qqline(resid(anova_result))

# 6.2. Extract the residuals and Homogeneity of Variances
aov_residuals <- residuals(object = anova_result )
aov_residuals <- resid(anova_result)
nortest::lillie.test(aov_residuals)
plot(anova_result, 1)
car::leveneTest(dN ~ any, data = data) #p-value is not less than the significance level of 0.05. This means that there is no evidence to suggest that the variance across groups is statistically significantly different. Therefore, we can assume the homogeneity of variances in the different treatment groups.
