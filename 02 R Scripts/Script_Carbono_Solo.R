### C NO SOLO: REGRESSÃO MISTO - SITE ALEATÓRIO ###

## DISTRIBUIÇÃO 
library(tidyverse)
library(writexl)
library(readxl)
setwd("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/distribuicao")
dadosmisto <- read_excel("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/dadosmisto.xlsx",row.names(1))


png("histograma_c_solo.png", width = 800, height = 600)
hist_biomass <- hist(dadosmisto$yeojohnson_c_soloid, main="Carbono no Solo", xlab="Valores", col="blue", border="black")
dev.off()
shapiro.test(dadosmisto$c_soloid) # p-value = 0.1823

dadosmisto <- dadosmisto %>% 
  mutate(log.c_soloid = scale(c_soloid))

png("histograma_c_solo.png", width = 800, height = 600)
hist_biomass <- hist(dadosmisto$log.c_soloid, main="Carbono no Solo", xlab="Valores", col="blue", border="black")
dev.off()
shapiro.test(dadosmisto$log.c_soloid) # p-value = 0.1823

str(dadosmisto$c_soloid)


dadosmisto <- dadosmisto %>%
  mutate(sqrt_c_soloid = sqrt(c_soloid))
shapiro.test(dadosmisto$sqrt_c_soloid)

library(MASS)
bc <- boxcox(lm(c_soloid ~ 1, data = dadosmisto), lambda = seq(-2, 2, 0.1))
lambda_opt <- bc$x[which.max(bc$y)] # Melhor lambda
dadosmisto <- dadosmisto %>%
  mutate(boxcox_c_soloid = ((c_soloid^lambda_opt - 1) / lambda_opt))
shapiro.test(dadosmisto$boxcox_c_soloid)

install.packages("bestNormalize")
library(bestNormalize)

yeo_johnson <- bestNormalize(dadosmisto$c_soloid, method = "yeojohnson")
dadosmisto$yeojohnson_c_soloid <- predict(yeo_johnson)

shapiro.test(dadosmisto$yeojohnson_c_soloid)

gamma_test(dadosmisto$c_soloid)

library(writexl)
library(readxl)
library(lmerTest) #MODELO MISTO
library(MuMIn) #AICc
setwd("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística")
dadosmisto <- read_excel("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/dadosmisto.xlsx")

# OFICIAL ### PRODUTIVIDADE ### site as random variable

# Modelo nulo sem os preditores fixos
modelo_nulo <- lmer(c_soloid ~ (1 | site), data = dadosmisto, REML = FALSE)
summary(modelo_nulo)
AICc(modelo_nulo)


misto.gini <- lmer(c_soloid ~ gini + (1 | site), data = dadosmisto, REML = FALSE)
summary(misto.gini) 
anova(modelo_nulo, misto.gini) # p-value = 0.00156 **
r.squaredGLMM(misto.gini) # 0.9493683
AICc(misto.gini) # -862.3059

misto.SR <- lmer(yeojohnson_c_soloid ~ SR + (1 | site), data = dadosmisto, REML = FALSE)
summary(misto.SR) 
anova(modelo_nulo, misto.SR) # p-value = 0.000233 ***
r.squaredGLMM(misto.SR) # 0.9490071
AICc(misto.SR) # -865.9235

misto.pse <- lmer(yeojohnson_c_soloid ~ PSE + (1 | site), data = dadosmisto, REML = FALSE)
summary(misto.pse) 
anova(modelo_nulo, misto.pse) # p-value = 0.9191
r.squaredGLMM(misto.pse) # 0.9389322
AICc(misto.pse) # -852.2665

misto.PD <- lmer(yeojohnson_c_soloid ~ PD + (1 | site), data = dadosmisto, REML = FALSE)
summary(misto.PD) # p-value = 7.81e-05 ***
r.squaredGLMM(misto.PD) # 0.9509005
AICc(misto.PD) # -867.9986

misto.sesMPD <- lmer(yeojohnson_c_soloid ~ SESMPD + (1 | site), data = dadosmisto, REML = FALSE)
summary(misto.sesMPD) # p-value = 0.02616 * 
r.squaredGLMM(misto.sesMPD) # 0.9445938
AICc(misto.sesMPD) # -857.2279

misto.sesMNTD <- lmer(yeojohnson_c_soloid ~ SESMNTD + (1 | site), data = dadosmisto, REML = FALSE)
summary(misto.sesMNTD) # p-value = 0.17537 
r.squaredGLMM(misto.sesMNTD) # 0.9398334
AICc(misto.sesMNTD) # -854.108

misto.PSV <- lmer(yeojohnson_c_soloid ~ PSV + (1 | site), data = dadosmisto, REML = FALSE)
summary(misto.PSV) # p-value = 0.001494 ** 
r.squaredGLMM(misto.PSV) # 0.9467189
AICc(misto.PSV) # -862.4304

misto.PSC <- lmer(yeojohnson_c_soloid ~ PSC + (1 | site), data = dadosmisto, REML = FALSE)
summary(misto.PSC) # p-value = 0.37578
r.squaredGLMM(misto.PSC) # 0.9390914
AICc(misto.PSC) # -853.047

misto.PPT <- lmer(yeojohnson_c_soloid ~ ppt + (1 | site), data = dadosmisto, REML = FALSE)
summary(misto.PPT) # p-value = 0.1780  
r.squaredGLMM(misto.PPT) # 0.9637944
AICc(misto.PPT) # -853.7271

misto.tmax <- lmer(yeojohnson_c_soloid ~ tmax + (1 | site), data = dadosmisto, REML = FALSE)
summary(misto.tmax) # p-value = 0.000261 ***  
r.squaredGLMM(misto.tmax) # 0.9747308
AICc(misto.tmax) # -864.2347

misto.tmin <- lmer(yeojohnson_c_soloid ~ tmin + (1 | site), data = dadosmisto, REML = FALSE)
summary(misto.tmin) # p-value = 0.0136 *  
r.squaredGLMM(misto.tmin) # 0.9687429
AICc(misto.tmin) # -857.5978

misto.pet <- lmer(yeojohnson_c_soloid ~ pet + (1 | site), data = dadosmisto, REML = FALSE)
summary(misto.pet) # p-value = 0.215  
r.squaredGLMM(misto.pet) # 0.9258293
AICc(misto.pet) # -853.4526

misto.vpd <- lmer(yeojohnson_c_soloid ~ vpd + (1 | site), data = dadosmisto, REML = FALSE)
summary(misto.vpd) # p-value = 0.00023 ***  
r.squaredGLMM(misto.vpd) # 0.978341
AICc(misto.vpd) # -864.0184

misto.alt <- lmer(yeojohnson_c_soloid ~ altitude + (1 | site), data = dadosmisto, REML = FALSE)
summary(misto.alt) # p-value = 2.98e-07 *** 
r.squaredGLMM(misto.alt) # 0.9903809
AICc(misto.alt) # -874.2442

misto.dec <- lmer(yeojohnson_c_soloid ~ declividade + (1 | site), data = dadosmisto, REML = FALSE)
summary(misto.dec) # p-value = 0.10104   
r.squaredGLMM(misto.dec) # 0.9294974
AICc(misto.dec) # -854.875

misto.simp <- lmer(yeojohnson_c_soloid ~ unbias.simp + (1 | site), data = dadosmisto, REML = FALSE)
summary(misto.simp) # p-value = 0.002937 **    
r.squaredGLMM(misto.simp) # 0.9451208
AICc(misto.simp) # -861.187

misto.shan <- lmer(yeojohnson_c_soloid ~ shannon + (1 | site), data = dadosmisto, REML = FALSE)
summary(misto.shan) # p-value = 0.000428 ***  
r.squaredGLMM(misto.shan) # 0.9478076
AICc(misto.shan) # -864.7851

misto.sr_inicial <- lmer(yeojohnson_c_soloid ~ riq_inicial + (1 | site), data = dadosmisto, REML = FALSE)
summary(misto.sr_inicial) # p-value = 0.19498  
r.squaredGLMM(misto.sr_inicial) # 0.9389875
AICc(misto.sr_inicial) # -853.9497

misto.pcps1 <- lmer(yeojohnson_c_soloid ~ pcps1 + (1 | site), data = dadosmisto, REML = FALSE)
summary(misto.pcps1) # p-value = 0.65090  

misto.pcps2 <- lmer(yeojohnson_c_soloid ~ pcps2 + (1 | site), data = dadosmisto, REML = FALSE)
summary(misto.pcps2) # p-value = 0.53232     

misto.PC1nutri <- lmer(yeojohnson_c_soloid ~ PC1nutri + (1 | site), data = dadosmisto, REML = FALSE)
summary(misto.PC1nutri) # p-value = 0.03948 *   
r.squaredGLMM(misto.PC1nutri) # 0.938592
AICc(misto.PC1nutri) # -856.5254

misto.PC2nutri <- lmer(yeojohnson_c_soloid ~ PC2nutri + (1 | site), data = dadosmisto, REML = FALSE)
summary(misto.PC2nutri) # p-value = 7.76e-06 ***   
r.squaredGLMM(misto.PC2nutri) # 0.9440901
AICc(misto.PC2nutri) # -872.5287

misto.valor_s <- lmer(yeojohnson_c_soloid ~ valor_s + (1 | site), data = dadosmisto, REML = FALSE)
summary(misto.valor_s) # p-value = 0.89516

misto.valor_t <- lmer(yeojohnson_c_soloid ~ valor_t + (1 | site), data = dadosmisto, REML = FALSE)
summary(misto.valor_t) # p-value = 6.06e-12 ***  
r.squaredGLMM(misto.valor_t) # 0.9662621
AICc(misto.valor_t) # -900.31

misto.ph <- lmer(yeojohnson_c_soloid ~ ph + (1 | site), data = dadosmisto, REML = FALSE)
summary(misto.ph) # p-value = 6.83e-05 ***  
r.squaredGLMM(misto.ph) # 0.9463835
AICc(misto.ph) # -868.2701

misto.faba <- lmer(yeojohnson_c_soloid ~ faba + (1 | site), data = dadosmisto, REML = FALSE)
summary(misto.faba) # p-value = 0.68781

misto.PC1 <- lmer(yeojohnson_c_soloid ~ PC1 + (1 | site), data = dadosmisto, REML = FALSE)
summary(misto.PC1) # p-value = 2.54e-13 *** 
r.squaredGLMM(misto.PC1) # 0.9671705
AICc(misto.PC1) # -906.7315

misto.PC2 <- lmer(yeojohnson_c_soloid ~ PC2 + (1 | site), data = dadosmisto, REML = FALSE)
summary(misto.PC2) # p-value = 1.07e-05 *** 
r.squaredGLMM(misto.PC2) # 0.9539789
AICc(misto.PC2) # -871.7946

misto.silte <- lmer(yeojohnson_c_soloid ~ silte_g.kg + (1 | site), data = dadosmisto, REML = FALSE)
summary(misto.silte) # p-value = 0.000162 ***
r.squaredGLMM(misto.silte) # 0.9499691
AICc(misto.silte) # -866.6022

misto.mcwd <- lmer(yeojohnson_c_soloid ~ mcwd + (1 | site), data = dadosmisto, REML = FALSE) # Some predictor variables are on very different scales: consider rescaling
summary(misto.mcwd) # p-value = 0.030394 * 
r.squaredGLMM(misto.mcwd) # 0.9555121
AICc(misto.mcwd) # -857.4783

#### MODELO GAMMA ####

library(fitdistrplus)
library(goftest)

# Ajustar uma distribuição Gama aos dados
fit_gamma <- fitdist(dadosmisto$c_soloid, "gamma")

# Verificar o ajuste visualmente
plot(fit_gamma)

# Teste de Kolmogorov-Smirnov
ks.test(dadosmisto$c_soloid, "pgamma", shape = fit_gamma$estimate["shape"], rate = fit_gamma$estimate["rate"])

# Teste de Anderson-Darling
ad.test(dadosmisto$c_soloid, "pgamma", shape = fit_gamma$estimate["shape"], rate = fit_gamma$estimate["rate"])

library(lme4)
modelo_gama <- glmer(c_soloid ~ gini + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_gama)

modelo_gama <- glmer(c_soloid ~ SR + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_gama)

modelo_gama <- glmer(c_soloid ~ PD + (1 | site), data = dadosmisto, family = Gamma(link = "log"))
summary(modelo_gama)

library(lme4)
modelo_lognormal <- lmer(log(c_soloid) ~ PD + (1 | site), data = dadosmisto)
summary(modelo_lognormal)

# Comparar AIC
AIC(modelo_gama, modelo_lognormal)



