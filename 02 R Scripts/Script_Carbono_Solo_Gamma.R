
install.packages("goftest")
library(fitdistrplus)
library(goftest)

# Ajustar uma distribuição Gama aos dados
fit_gamma <- fitdist(dadosmisto$c_soloid, "gamma")
fit_gamma <- fitdist(dadosmisto$n_sid, "gamma")
fit_gamma <- fitdist(dadosmisto$c.n_soloid, "gamma")

# Verificar o ajuste visualmente
plot(fit_gamma)

# Teste de Anderson-Darling
ad.test(dadosmisto$c_soloid, "pgamma", shape = fit_gamma$estimate["shape"], rate = fit_gamma$estimate["rate"])
ad.test(dadosmisto$n_sid, "pgamma", shape = fit_gamma$estimate["shape"], rate = fit_gamma$estimate["rate"])
ad.test(dadosmisto$c.n_soloid, "pgamma", shape = fit_gamma$estimate["shape"], rate = fit_gamma$estimate["rate"])

library(ggplot2)
ggplot(dadosmisto, aes(x = n_sid)) +
  geom_density(fill = "blue", alpha = 0.5) +
  stat_function(fun = dgamma, args = list(shape = fit_gamma$estimate["shape"], rate = fit_gamma$estimate["rate"]), color = "red", size = 1) +
  labs(title = "Ajuste da distribuição Gama", x = "n_sid", y = "Densidade")

ggplot(dadosmisto, aes(x = c.n_soloid)) +
  geom_density(fill = "blue", alpha = 0.5) +
  stat_function(fun = dgamma, args = list(shape = fit_gamma$estimate["shape"], rate = fit_gamma$estimate["rate"]), color = "red", size = 1) +
  labs(title = "Ajuste da distribuição Gama", x = "n_sid", y = "Densidade")



# C no solo

library(writexl)
library(readxl)
library(MuMIn) #AICc
library(lme4)

modelo_gini <- glmer(c_soloid ~ gini + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_gini) # p-value = 0.205    

modelo_sr <- glmer(c_soloid ~ SR + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_sr) # p-value = 0.356     

modelo_pd <- glmer(c_soloid ~ SESPD + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_pd) # p-value = 0.946 

modelo_sesmpd <- glmer(c_soloid ~ SESMPD + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_sesmpd) # p-value = 0.902 

modelo_sesmntd <- glmer(c_soloid ~ SESMNTD + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_sesmntd) # p-value = 0.999

modelo_PSV <- glmer(c_soloid ~ PSV + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_PSV) # p-value = 0.779 

modelo_PSC <- glmer(c_soloid ~ PSC + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_PSC) # p-value = 0.569

modelo_PSE <- glmer(c_soloid ~ PSE + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_PSE) # p-value = 0.679

modelo_ppt <- glmer(c_soloid ~ ppt + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_ppt) # p-value = 0.000486 *** AIC -244.5

modelo_tmax <- glmer(c_soloid ~ tmax + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_tmax) # p-value = 6.75e-08 *** AIC -256.0

modelo_tmin <- glmer(c_soloid ~ tmin + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_tmin) # p-value = 4.37e-06 *** AIC -250.4

modelo_pet <- glmer(c_soloid ~ pet + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_pet) # p-value = 0.617

modelo_vpd <- glmer(c_soloid ~ vpd + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_vpd) # p-value = 6.28e-08 *** AIC -255.9

modelo_alt <- glmer(c_soloid ~ alt_scale + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_alt) # p-value = 7.15e-05 *** AIC -246.7

modelo_dec <- glmer(c_soloid ~ declividade + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_dec) # p-value = 0.208

modelo_simp <- glmer(c_soloid ~ unbias.simp + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_simp) # p-value = 0.445

modelo_shan <- glmer(c_soloid ~ shannon + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_shan) # p-value = 0.371    

dadosmisto$mcwd_scale <- scale(dadosmisto$mcwd)
modelo_mcwd <- glmer(c_soloid ~ mcwd_scale + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_mcwd) # p-value = 0.846

modelo_pcps1 <- glmer(c_soloid ~ pcps1 + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_pcps1) # p-value = 0.104    

modelo_pcps2 <- glmer(c_soloid ~ pcps2 + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_pcps2) # p-value = 0.748     

dadosmisto$riq_inicial_scaled <- scale(dadosmisto$riq_inicial)
modelo_sri_scaled <- glmer(c_soloid ~ riq_inicial_scaled + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_sri_scaled) # p-value = 0.0814 . AIC -235.8

modelo_PC1nutri <- glmer(c_soloid ~ PC1nutri + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_PC1nutri) # p-value = 0.0023 ** AIC -241.0
r.squaredGLMM(modelo_PC1nutri)

modelo_PC2nutri <- glmer(c_soloid ~ PC2nutri + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_PC2nutri) # p-value = 0.301

modelo_valor_s <- glmer(c_soloid ~ valor_s + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_valor_s) # p-value = 0.707    

modelo_valor_t <- glmer(c_soloid ~ valor_t + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_valor_t) # p-value = 2.51e-09 *** AIC -261.9

modelo_ph <- glmer(c_soloid ~ ph + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_ph) # p-value = 0.8383  

modelo_faba <- glmer(c_soloid ~ faba + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_faba) # p-value = 0.117

dadosmisto$silte_scale <- scale(dadosmisto$silte_g.kg)
modelo_silte <- glmer(c_soloid ~ silte_scale + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_silte) # p-value = 0.015 * AIC -238.5

modelo_PC1 <- glmer(c_soloid ~ PC1 + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_PC1) # p-value = 2.31e-09 *** AIC -261.9
r.squaredGLMM(modelo_PC1)

modelo_PC2 <- glmer(c_soloid ~ PC2 + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_PC2) # p-value = 0.0904 . AIC -235.7

modelo_geral <- glmer(c_soloid ~ valor_t + riq_inicial_scaled + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_geral) # AIC -268.9 

modelo_geral <- glmer(c_soloid ~ valor_t + SR + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_geral) # AIC -263.9

modelo_geral <- glmer(c_soloid ~ PC1 + SR + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_geral) # AIC -266.8

modelo_geral <- glmer(c_soloid ~ valor_t + riq_inicial_scaled + SR + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_geral) # AIC -270.7 
r.squaredGLMM(modelo_geral)

modelo_geral <- glmer(c_soloid ~ PC1 + SR + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_geral) # AIC -266.8 

modelo_interacao <- glmer(c_soloid ~ valor_t + riq_inicial_scaled + SR + SR:valor_t + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_interacao)

cor(dadosmisto[, c("valor_t", "riq_inicial_scaled", "SR", "PC2")])

# N no solo

library(lme4)
modelo_gini <- glmer(n_sid ~ gini + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_gini) # p-value = 0.205    

modelo_SR <- glmer(n_sid ~ SR + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_SR) # p-value = 0.302    

modelo_PD <- glmer(n_sid ~ SESPD + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_PD) # p-value = 0.636     

modelo_MPD <- glmer(n_sid ~ SESMPD + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_MPD) # p-value = 0.529

modelo_MNTD <- glmer(n_sid ~ SESMNTD + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_MNTD) # p-value = 0.603

modelo_PSV <- glmer(n_sid ~ PSV + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_PSV) # p-value = 0.589

modelo_PSC <- glmer(n_sid ~ PSC + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_PSC) # p-value = 0.265 

modelo_PSE <- glmer(n_sid ~ PSE + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_PSE) # p-value = 0.31

modelo_pcps1 <- glmer(n_sid ~ pcps1 + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_pcps1) # p-value = 0.0266 *

modelo_pcps2 <- glmer(n_sid ~ pcps2 + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_pcps2) # p-value = 0.935 

dadosmisto$shan_scaled <- scale(dadosmisto$shannon)
modelo_shannon <- glmer(n_sid ~ shan_scaled + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_shannon) # p-value = 0.00219 ** 

modelo_simp <- glmer(n_sid ~ unbias.simp + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_simp) # p-value = 0.217

modelo_PC1nutri <- glmer(n_sid ~ PC1nutri + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_PC1nutri) # p-value = 0.172  

modelo_faba <- glmer(n_sid ~ faba + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_faba) # p-value = 0.0514 . 

dadosmisto$silte_scaled <- scale(dadosmisto$silte_g.kg)
modelo_silte <- glmer(n_sid ~ silte_scaled + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_silte) # p-value = 9.06e-07 ***

modelo_PC1 <- glmer(n_sid ~ PC1 + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_PC1) # p-value = 3.66e-09 ***

modelo_PC2 <- glmer(n_sid ~ PC2 + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_PC2) # p-value = 0.667    

dadosmisto$ppt_scale <- scale(dadosmisto$ppt)
modelo_PPT <- glmer(n_sid ~ ppt_scale + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_PPT) # p-value = 0.734

modelo_tmax <- glmer(n_sid ~ tmax + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_tmax) # p-value = 2.93e-05 ***

modelo_tmin <- glmer(n_sid ~ tmin + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_tmin) # p-value =  0.0162 * 

dadosmisto$pet_scale <- scale(dadosmisto$pet)
modelo_pet <- glmer(n_sid ~ pet_scale + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_pet) # p-value =  0.142

modelo_vpd <- glmer(n_sid ~ vpd + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_vpd) # p-value =  2.07e-05 ***

dadosmisto$alt_scale <- scale(dadosmisto$altitude)
modelo_alt <- glmer(n_sid ~ alt_scale + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_alt) # p-value = 4.78e-08 ***

modelo_dec <- glmer(n_sid ~ declividade + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_dec) # p-value = 0.319

modelo_mcwd <- glmer(n_sid ~ mcwd_scale + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_mcwd) # p-value = 0.699

# C/N no solo

library(lme4)
modelo_gini <- glmer(c.n_soloid ~ gini + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_gini) # p-value = 0.18876       

modelo_SR <- glmer(c.n_soloid ~ SR + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_SR) # p-value = 0.00502 **  AIC 86.6   

modelo_PD <- glmer(c.n_soloid ~ SESPD + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_PD) # p-value = 0.468200         

modelo_MPD <- glmer(c.n_soloid ~ SESMPD + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_MPD) # p-value = 0.599257    

modelo_MNTD <- glmer(c.n_soloid ~ SESMNTD + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_MNTD) # p-value = 0.480426    

modelo_PSV <- glmer(c.n_soloid ~ PSV + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_PSV) # p-value = 0.3800  

modelo_PSC <- glmer(c.n_soloid ~ PSC + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_PSC) # p-value = 0.954735     

modelo_PSE <- glmer(c.n_soloid ~ PSE + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_PSE) # p-value = 0.928718    

modelo_pcps1 <- glmer(c.n_soloid ~ pcps1 + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_pcps1) # p-value = 0.308419

modelo_pcps2 <- glmer(c.n_soloid ~ pcps2 + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_pcps2) # p-value = 0.757667 

dadosmisto$shan_scaled <- scale(dadosmisto$shannon)
modelo_shannon <- glmer(c.n_soloid ~ shan_scaled + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_shannon) # p-value = 0.015125 * AIC 88.4 

modelo_simp <- glmer(c.n_soloid ~ unbias.simp + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_simp) # p-value = 0.0882 . AIC 91.2    

modelo_PC1nutri <- glmer(c.n_soloid ~ PC1nutri + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_PC1nutri) # p-value = 0.0857 .  AIC 91.2    

modelo_faba <- glmer(c.n_soloid ~ faba + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_faba) # p-value = 0.325937     

dadosmisto$silte_scaled <- scale(dadosmisto$silte_g.kg)
modelo_silte <- glmer(c.n_soloid ~ silte_scaled + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_silte) # p-value = 0.000426 ***

modelo_PC1 <- glmer(c.n_soloid ~ PC1 + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_PC1) # p-value = 0.636140 

modelo_PC2 <- glmer(c.n_soloid ~ PC2 + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_PC2) # p-value = 0.50013        

dadosmisto$ppt_scale <- scale(dadosmisto$ppt)
modelo_PPT <- glmer(c.n_soloid ~ ppt_scale + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_PPT) # p-value = 0.031663 * AIC 86.8     

modelo_tmax <- glmer(c.n_soloid ~ tmax + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_tmax) # p-value = 0.188

modelo_tmin <- glmer(c.n_soloid ~ tmin + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_tmin) # p-value = 0.026 * 

dadosmisto$pet_scale <- scale(dadosmisto$pet)
modelo_pet <- glmer(c.n_soloid ~ pet_scale + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_pet) # p-value =  0.095239 . AIC 90.7

modelo_vpd <- glmer(c.n_soloid ~ vpd + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_vpd) # p-value = 0.252

dadosmisto$alt_scale <- scale(dadosmisto$altitude)
modelo_alt <- glmer(c.n_soloid ~ alt_scale + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_alt) # p-value = 0.45148   

modelo_dec <- glmer(c.n_soloid ~ declividade + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_dec) # p-value = 1.39e-05 *** AIC 92.8

modelo_mcwd <- glmer(c.n_soloid ~ mcwd_scale + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_mcwd) # p-value = 0.5360  

################
### GRÁFICOS ###
################

library(tidyverse)
library(ggplot2)
setwd("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/graficos")

## valor-t 
dadosmisto %>%
  ggplot(aes(valor_t,c_soloid))+
  geom_point(size=3,alpha=0.5)+ 
  geom_smooth(method = lm,se=T,colour="red")+
  labs(x="t-value",
       y="Soil Carbon (CHN%)",
       title = "Relationship between t-value and Soil Carbon")+
  annotate("text",x=11,y=0.17,label="R²c = 0.76")+
  theme_replace()

ggsave("csoloid_valor.t.jpeg",width = 15,height = 10,units = "cm")

## PC1 
dadosmisto %>%
  ggplot(aes(PC1,c_soloid))+
  geom_point(size=3,alpha=0.5)+ 
  geom_smooth(method = lm,se=T,colour="purple")+
  labs(x="PCA - granulometry - axis 1",
       y="Soil Carbon (CHN%)",
       title = "Relationship between Granulometry and Soil Carbon")+
  annotate("text",x=3,y=0.17,label="R²c = 0.67")+
  theme_replace()

ggsave("csoloid_PC1.jpeg",width = 15,height = 10,units = "cm")

## valor-t 
dadosmisto %>%
  ggplot(aes(PC1nutri,c_soloid))+
  geom_point(size=3,alpha=0.5)+ 
  geom_smooth(method = lm,se=T,colour="yellow")+
  labs(x="PCA soil nutrients - axis 1",
       y="Soil Carbon (CHN%)",
       title = "Relationship between soil nutrients and Soil Carbon")+
  annotate("text",x=1.5,y=0.17,label="R²c = 0.53")+
  theme_replace()

ggsave("csoloid_PC1nutri.jpeg",width = 15,height = 10,units = "cm")

## riqueza
dadosmisto %>%
  ggplot(aes(SR,c_soloid))+
  geom_point(size=3,alpha=0.5)+ 
  geom_smooth(method = lm,se=T,colour="yellow")+
  labs(x="PCA soil nutrients - axis 1",
       y="Soil Carbon (CHN%)",
       title = "Relationship between soil nutrients and Soil Carbon")+
  annotate("text",x=1.5,y=0.17,label="R²c = 0.53")+
  theme_replace()
