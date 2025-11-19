# Importar o arquivo de dados em formato .csv

regressao.produt <- read.csv("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Dados/dadosreg1.csv", header = T, sep = ";")
str(biomassa)

library(tidyverse)
library(ggplot2)

### log PRODUT ###

regressao.produt <- regressao.produt %>%
  mutate(logprodut_z = log(produt_gm2ano_z))



### PRODUTIVIDADE Zanini G/M2/ANO OFICIAL SR ###
regressao.produt %>%
  ggplot(aes(SR,logprodut_z))+
  geom_point(size=3,alpha=0.5)+ 
  geom_smooth(method = lm,se=T,colour="orange")+
  labs(x="SR",
       y="log Produtividade (g/m²/ano)",
       title = "Produtividade em função da Riqueza de espécies")+
  theme_replace()+
  annotate("text",x=15,y=5.5,label="R² = 0.207 
p-valor = 0.0007")

ggsave("ofc_logprodutzanini_sr1.png",width = 15,height = 10,units = "cm")

logprodutz.SR <- lm(logprodut_z~SR,data = regressao.produt)
summary(logprodutz.SR)


### PRODUTIVIDADE Zanini G/M2/ANO OFICIAL PD ###
regressao.produt %>%
  ggplot(aes(PD,logprodut_z))+
  geom_point(size=3,alpha=0.5)+ 
  geom_smooth(method = lm,se=T,colour="orange")+
  labs(x="PD",
       y="log Produtividade (g/m²/ano)",
       title = "Produtividade em função do PD Faith")+
  theme_replace()+
  annotate("text",x=1250,y=5.5,label="R² = 0.144 
p-valor = 0.005")

ggsave("ofc_logprodutzanini_pd.png",width = 15,height = 10,units = "cm")

logprodutz.PD <- lm(logprodut_z~PD,data = regressao.produt)
summary(logprodutz.PD)

### PRODUTIVIDADE Zanini G/M2/ANO OFICIAL ses.PD ###
regressao.produt %>%
  ggplot(aes(SES.PD,logprodut_z))+
  geom_point(size=3,alpha=0.5)+ 
  geom_smooth(method = lm,se=T,colour="orange")+
  labs(x="SES.PD",
       y="log Produtividade (g/m²/ano)",
       title = "Produtividade em função de SES.PD Faith")+
  theme_replace()+
  annotate("text",x=0.5,y=5.5,label="R² = 0.008 
p-valor = 0.534")

ggsave("ofc_logprodutzanini_sespd1.png",width = 15,height = 10,units = "cm")

logprodutz.SESPD <- lm(logprodut_z~SES.PD,data = regressao.produt)
summary(logprodutz.SESPD)

### PRODUTIVIDADE Zanini G/M2/ANO OFICIAL MPD ###
regressao.produt %>%
  ggplot(aes(MPD,logprodut_z))+
  geom_point(size=3,alpha=0.5)+ 
  geom_smooth(method = lm,se=T,colour="orange")+
  labs(x="MPD",
       y="log Produtividade (g/m²/ano)",
       title = "Produtividade em função do MPD",
       subtitle = "Segundo equação alométrica de Zanini et al., 2022")+
  theme_replace()+
  annotate("text",x=0.65,y=5.5,label="R² = 0.002 
p-valor = 0.7154")

ggsave("ofc_logprodutzanini_MPD.png",width = 15,height = 10,units = "cm")

logprodutz.MPD <- lm(logprodut_z~MPD,data = regressao.produt)
summary(logprodutz.MPD)

### PRODUTIVIDADE Zanini G/M2/ANO OFICIAL ses.MPD ###
regressao.produt %>%
  ggplot(aes(SES.MPD,logprodut_z))+
  geom_point(size=3,alpha=0.5)+ 
  geom_smooth(method = lm,se=T,colour="orange")+
  labs(x="SES.MPD",
       y="log Produtividade (g/m²/ano)",
       title = "Produtividade em função do MPD padronizado pela riqueza",
       subtitle = "Segundo equação alométrica de Zanini et al., 2022")+
  theme_replace()+
  annotate("text",x=-4.5,y=5.5,label="R² = 0.0178 
p-valor = 0.345")

ggsave("ofc_logprodutzanini_sesmpd.png",width = 15,height = 10,units = "cm")

logprodutz.SESMPD <- lm(logprodut_z~SES.MPD,data = regressao.produt)
summary(logprodutz.SESMPD)

### PRODUTIVIDADE Zanini G/M2/ANO OFICIAL ses.MNTD ###
regressao.produt %>%
  ggplot(aes(SES.MNTD,logprodut_z))+
  geom_point(size=3,alpha=0.5)+ 
  geom_smooth(method = lm,se=T,colour="orange")+
  labs(x="SES.MNTD",
       y="log Produtividade (g/m²/ano)",
       title = "Produtividade em função do MNTD padronizado pela riqueza",
       subtitle = "Segundo equação alométrica de Zanini et al., 2022")+
  theme_replace()+
  annotate("text",x=1,y=5.5,label="R² = 0.0001 
p-valor = 0.938")

ggsave("ofc_logprodutzanini_sesmntd.png",width = 15,height = 10,units = "cm")

logprodutz.SESMNTD <- lm(logprodut_z~SES.MNTD,data = regressao.produt)
summary(logprodutz.SESMNTD)


### PRODUTIVIDADE Zanini G/M2/ANO OFICIAL MNTD ###
regressao.produt %>%
  ggplot(aes(MNTD,logprodut_z))+
  geom_point(size=3,alpha=0.5)+ 
  geom_smooth(method = lm,se=T,colour="orange")+
  labs(x="MNTD",
       y="log Produtividade (g/m²/ano)",
       title = "Produtividade em função do MNTD",
       subtitle = "Segundo equação alométrica de Zanini et al., 2022")+
  theme_replace()+
  annotate("text",x=0.75,y=5.5,label="R² = 0.002 
p-valor = 0.7307")

ggsave("ofc_logprodutzanini_mntd.png",width = 15,height = 10,units = "cm")

logprodutz.MNTD <- lm(logprodut_z~MNTD,data = regressao.produt)
summary(logprodutz.MNTD)
#####################################################



#######
# PSV #
#######

library(tidyverse)

### ZANINI E PSV ###
regressao.produt %>%
  ggplot(aes(PSV,logprodut_z))+
  geom_point(size=3,alpha=0.5)+ 
  geom_smooth(method = lm,se=T,colour="orange")+
  labs(x="PSV",
       y="log Produtividade (g/m²/ano)",
       title = "Produtividade em função do PSV",
       subtitle = "Segundo equação alométrica de Zanini et al., 2022")+
  theme_replace()+
  annotate("text",x=0.6,y=5.5,label="R² = 0.00039 
p-valor = 0.8888")

ggsave("ofc_logprodutzanini_PSV.png",width = 15,height = 10,units = "cm")

logprodutz.PSV <- lm(logprodut_z~PSV,data = regressao.produt)
summary(logprodutz.PSV)

### ZANINI E PSC ###
regressao.produt %>%
  ggplot(aes(PSC,logprodut_z))+
  geom_point(size=3,alpha=0.5)+ 
  geom_smooth(method = lm,se=T,colour="orange")+
  labs(x="PSC",
       y="log Produtividade (g/m²/ano)",
       title = "Produtividade em função do PSC",
       subtitle = "Segundo equação alométrica de Zanini et al., 2022")+
  theme_replace()+
  annotate("text",x=0.6,y=5.5,label="R² = 0.025 
p-valor = 0.2623")

ggsave("ofc_logprodutzanini_PSC.png",width = 15,height = 10,units = "cm")

logprodutz.PSC <- lm(logprodut_z~PSC,data = regressao.produt)
summary(logprodutz.PSC)

### ZANINI E PSR ###
regressao.produt %>%
  ggplot(aes(PSR,logprodut_z))+
  geom_point(size=3,alpha=0.5)+ 
  geom_smooth(method = lm,se=T,colour="orange")+
  labs(x="PSR",
       y="log Produtividade (g/m²/ano)",
       title = "Produtividade em função da Phylogenetic Species Richness")+
  theme_replace()+
  annotate("text",x=11.2,y=5.5,label="R² = 0.172 
p-valor = 0.002")

ggsave("ofc_logprodutzanini_sesPSR1.png",width = 15,height = 10,units = "cm")

logprodutz.PSR <- lm(logprodut_z~PSR,data = regressao.produt)
summary(logprodutz.PSR)
