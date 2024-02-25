#Fer plot Beta- coefficients
require(tidyverse)
set.seed(123)
dat <- data.frame(x = rnorm(100), z = rnorm(100), y1 = rnorm(100), y2 = rnorm(100), y3 = rnorm(100))
mod1 <- m1
mod2 <- m2


## Create data frame of model coefficients and standard errors
# Function to extract what we need
ce = function(model.obj) {
  extract = summary(get(model.obj))$coefficients[ ,1:2]
  return(data.frame(extract, vars=row.names(extract), model=model.obj))
}

# Run function on the three models and bind into single data frame
coefs = do.call(rbind, sapply(paste0("mod",1:2), ce, simplify=FALSE))

names(coefs)[2] = "se" 

# Using $ notation
coefs$vars

Indexs = c( "(Intercept)","poly(julian_day, 2)1", "poly(julian_day, 2)2", "NAO index","AMO index" ,"(Intercept)", "poly(julian_day, 2)1", "poly(julian_day, 2)2", "AMO index")
coefs$Indexs = Indexs

coefs<-coefs[-1,]
coefs<-coefs[-1,]
coefs<-coefs[-1,]
coefs<-coefs[-3,]
coefs<-coefs[-3,]
coefs<-coefs[-3,]
# Faceted coefficient plot
coefs_mod1 <- filter(coefs, model== "mod1")
coefs_mod2 <- filter(coefs, model== "mod2")
pd <- position_dodge(width = 0.9)
mypalette<-brewer.pal(9,"Pastel1")

coefs$Indexs <- factor(coefs$Indexs, levels = unique(coefs$Indexs))

## Reverse the order of Subclass_Name levels
coefs$Indexs <- factor(coefs$Indexs,
                                 levels=c("NAO index","AMO index"))

levels=c("NAO index","AMO index")

coefs$Indexs <- factor(coefs$Indexs, levels=rev(c("NAO index","AMO index")))

ggplot2::ggplot(coefs, aes(Indexs, Estimate, group=model)) + 
  geom_hline(yintercept=0, lty=2, lwd=1, colour="lightgrey") +
  geom_errorbar(aes(ymin=Estimate - se, ymax=Estimate + se), width=0.2, position = position_dodge(width=0.5)) +
  geom_point(aes(shape=model, fill=model), position = position_dodge(0.5),size=3, color='black') +
  coord_flip() +
  scale_x_discrete(breaks=levels,labels=c("Avg NAO index","Avg AMO index")) +
  scale_shape_manual(values = c(21:22),labels = c(expression(paste(delta^15, "N")), expression(paste(delta^13, "C")))) +
  scale_fill_manual(values=c('#FC8D59','#74A9CF'),labels = c(expression(paste(delta^15, "N")), expression(paste(delta^13, "C")))) +
  guides(colour=FALSE) +
  labs(x="", y="Standardized β coefficients") +
  theme_article(base_size=20) + theme(legend.title = element_blank(),legend.position = c(.9, .9),aspect.ratio=3/4)

#plot time windows
#Amb julian day, 60-0 months

#dN

a<-head(Bp_dN_sw_weather_m_rel[[1]]$Dataset, 20) ## Lin_func Mean NAO, close window-open window (34-31) (beta coefficient= -0.3619122)
a <- a %>% filter(deltaAICc < (a$deltaAICc[1])+2)
nrow_a <- nrow(a)
a$Index <- c(rep("NAO", nrow_a))
b<-head(Bp_dN_sw_weather_m_rel[[2]]$Dataset, 20) ## Lin_func Mean AMO, close window-open window (52-35) (beta coefficient= -4.377916)
b <- b %>% filter(deltaAICc < (a$deltaAICc[1])+2)
nrow_b <- nrow(b)
b$Index <- c(rep("AMO", nrow_b))

data <- rbind(a,b)

levels=c("NAO","AMO")

data$Index <- factor(data$Index, levels=rev(c("NAO","AMO")))

m1 <- ggplot(data) +
  geom_vline(aes(xintercept = data$Index), color = "gray", linetype = "dotted") +
  geom_linerange(aes(x = data$Index, ymin = WindowOpen, ymax = WindowClose), color = '#74A9CF', size = 20, alpha = 0.9) +  # Afegeix les línies de separació
  geom_hline(yintercept = seq(1,48,1), lty = 1, lwd = 0.1, colour = "white") +
  geom_hline(yintercept = 0, lty = 1, lwd = 1, colour = "black") +
  coord_flip() +
  labs(x = "", y = "Months from sample", title = "") +
  theme_article(base_size = 20) + ylim(36,0)+
  scale_x_discrete(breaks=levels,labels=c("Avg NAO index","Avg AMO index")) + theme(aspect.ratio=3/4)

#dC

a1<-head(Bp_dC_sw_weather_m_rel[[2]]$Dataset, 20) ## Lin_func Mean NAO, close window-open window (34-31) (beta coefficient= -0.3619122)
a1 <- a1 %>% filter(deltaAICc < (a1$deltaAICc[1])+2)
nrow_a1 <- nrow(a1)
a1$Index <- c(rep("AMO", nrow_a1))


data1 <- rbind(a1)

levels1=c("AMO")

data1$Index <- factor(data1$Index, levels1=c(levels1))

m2 <- ggplot(data1) +
  geom_vline(aes(xintercept = data1$Index), color = "gray", linetype = "dotted") +
  geom_linerange(aes(x = data1$Index, ymin = WindowOpen, ymax = WindowClose), color = '#FC8D59', size = 20, alpha =0.9) +  # Afegeix les línies de separació
  geom_hline(yintercept = seq(1,48,1), lty = 1, lwd = 0.1, colour = "white") +
  geom_hline(yintercept = 0, lty = 1, lwd = 1, colour = "black") +
  coord_flip() +
  labs(x = "", y = "Months from sample", title = "") +
  theme_article(base_size = 20) + ylim(36,0)+
  scale_x_discrete(breaks=levels1,labels=c("Avg AMO index")) + theme(aspect.ratio=3/4)

m1+m2

m3 <- ggplot(data) +
  geom_vline(aes(xintercept = data$Index), color = "gray", linetype = "dotted") +
  geom_linerange(aes(x = data$Index, ymin = WindowOpen, ymax = WindowClose), color = '#74A9CF', size = 20, alpha = 0.9) +  # Afegeix les línies de separació
  geom_hline(yintercept = seq(1,48,1), lty = 1, lwd = 0.1, colour = "white") +
  geom_hline(yintercept = 0, lty = 1, lwd = 1, colour = "black") +
  geom_linerange(data=data1, aes(x = Index, ymin = WindowOpen, ymax = WindowClose), color = '#FC8D59', size = 20, alpha =0.9) +  # Afegeix les línies de separació
  coord_flip() +
  labs(x = "", y = "Months from sample", title = "") +
  theme_article(base_size = 20) + ylim(36,0)+
  scale_x_discrete(breaks=levels,labels=c("Avg NAO index","Avg AMO index")) + theme(aspect.ratio=3/4)
m3
