
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

# 2. Import data

data <- read_excel("~/Desktop/Doctorat/Analisis_isotops_barbes/Clima/Data/Krill_Islandia.xlsx", sheet = "Hoja2")

# 3. Boxplot for nitrogen
krill_N <- ggplot(data, aes(as.factor(any),dN)) + 
  geom_boxplot(fill="#74A9CF", alpha=0.5, width=0.6) + 
  ylim(4,9) +
  geom_jitter(color="black", size=1, alpha=0.9, shape=1, width = 0.1) +
  #stat_summary(fun.y=mean, geom="point", shape=20, size=2, color="red", fill="red") +
  xlab("Year") + 
  ylab(expression(paste(delta^{15}, "N (‰)"))) + 
  theme_article(base_size = 20) + theme(aspect.ratio = 3/4)

ggsave("/Users/marcruizisagales/Documents/GitHub/Climate-baleen-plate-isotopes/png/Figure_5_krill_Nitrogen_boxplot.png", krill_N, 
       device = png(width = 600, height = 400)
)

# 4. Boxplot for carbon
ggplot(data, aes(as.factor(any),dC)) + 
  geom_boxplot(fill="#74A9CF", alpha=0.5, width=0.6) +
  geom_jitter(color="black", size=1, alpha=0.9, shape=1, width = 0.1) +
  #stat_summary(fun.y=mean, geom="point", shape=20, size=2, color="red", fill="red") +
  xlab("Year") + 
  ylab(expression(paste(delta^{15}, "C (‰)"))) + 
  theme_article(base_size = 20) 


# 5. Perform ANOVA
data$any <- as.factor(data$any)
anova_result <- aov(dN ~ any, data = data)
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
