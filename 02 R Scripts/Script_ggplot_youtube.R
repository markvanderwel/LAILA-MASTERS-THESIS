### data visualisation

# data
# mapping (aesthetics): se vai estar no x ou y, a cor etc
# geometric representation (boxplot, dotplot,linegraph, barchat, histogram)
# statistics
# facet
# coordinate space
# labels
# theme

install.packages("tidyverse")
library(tidyverse)
library(ggplot2)

data()
?BOD

ggplot(data = BOD,mapping = aes(x=Time,
                                y=demand))+
  geom_point(size=5)+
  geom_line(colour="red")

# Outra forma de fazer, mais limpa
ggplot(BOD,aes(Time,demand))+
  geom_point(size=5)+
  geom_line(colour="red")

# Outro dataset

View(CO2)

# com o %>% o ggplot já sabe qual data tô usando
CO2 %>%
  ggplot(aes(conc,uptake,
             colour = Treatment))+
  geom_point(size=3,alpha=0.5)+ #alpha é transparência
  geom_smooth(method = lm,se=F)+ #linear model, se = F não tem intervalo de confi
  facet_wrap(~Type)+ # divide em dois gráficos diferentes
  labs(title = "Concentration of CO2")+
  theme_bw()



CO2 %>%
  ggplot(aes(Treatment,uptake))+
  geom_boxplot()+
  geom_point(alpha=0.5,
    aes(size = conc,
                 colour=Plant))+
  facet_wrap(~Type)+
  coord_flip()+
  theme_bw()+
  labs(title = "Chilled vs Non-chilled")

View(mpg)

mpg %>%
  filter(cty<25) %>%
  ggplot(aes(displ,cty))+
  geom_point(aes(colour = drv,
                 size = trans),
             alpha=0.5)+
  geom_smooth(method = lm)+
  facet_wrap(~year,nrow = 1)+
  labs(x="Engine size",
       y="MPH",
       title = "Fuel efficiency")+
  theme_bw()
  

names(mpg)

mpg %>%
  filter(cty<35) %>%
  ggplot(aes(displ,hwy))+
  geom_point(aes(colour = drv))+
  geom_smooth(method = lm,se=F)+
  labs(x="Engine size",
       y="MPG on the Highway",
       title = "Fuel efficiency")+
  theme_minimal()



