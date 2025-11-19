
dadosreg_s <- read_excel("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/dadosreg_solo.xlsx",row.names(1))

rm(medias_por_site)

library(dplyr)

medias_sites <- dadosreg_s %>%
  group_by(site) %>%  # Agrupa por site
  summarise(across(where(is.numeric), mean, na.rm = TRUE))  # Calcula a média apenas de colunas numéricas

print(medias_sites)

write_xlsx(medias_sites, "C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/dadosreg_solo1.xlsx")

library(writexl)
library(readxl)
setwd("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/distribuicao")
dadosreg_s <- read_excel("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/dadosreg_solo1.xlsx",row.names(1))


png("histograma_c_solo.png", width = 800, height = 600)
hist_biomass <- hist(dadosreg_s$c_solo, main="Carbono no Solo", xlab="Valores", col="blue", border="black")
dev.off()
shapiro.test(dadosreg_s$c_solo) # p-value = 0.1823

png("histograma_c_soloid.png", width = 800, height = 600)
hist_biomass <- hist(dadosreg_s$c_sid, main="Carbono no Solo - idade", xlab="Valores", col="purple", border="black")
dev.off()
shapiro.test(dadosreg_s$c_sid) # p-value = 0.5327

png("histograma_c_soloid.png", width = 800, height = 600)
hist_biomass <- hist(dadosmisto$c_soloid, main="Carbono no Solo - idade", xlab="Valores", col="purple", border="black")
dev.off()
shapiro.test(dadosmisto$c_soloid) # p-value = 0.5327

png("histograma_n_solo.png", width = 800, height = 600)
hist_biomass <- hist(dadosreg_s$n_solo, main="Nitrogênio no Solo", xlab="Valores", col="orange", border="black")
dev.off()
shapiro.test(dadosreg_s$n_solo) # p-value = 0.007878 (não é normal)

png("histograma_n_soloid.png", width = 800, height = 600)
hist_biomass <- hist(dadosreg_s$n_sid, main="Nitrogênio no Solo - idade", xlab="Valores", col="yellow", border="black")
dev.off()
shapiro.test(dadosreg_s$n_sid) # p-value = 0.111

png("histograma_c.n_solo.png", width = 800, height = 600)
hist_biomass <- hist(dadosreg_s$c.n_solo, main="C/N no Solo", xlab="Valores", col="brown", border="black")
dev.off()
shapiro.test(dadosreg_s$c.n_solo) # p-value = 0.08929

png("histograma_c.n_id.png", width = 800, height = 600)
hist_biomass <- hist(dadosreg_s$c.n_id, main="C/N no Solo - idade", xlab="Valores", col="blue", border="black")
dev.off()
shapiro.test(dadosreg_s$c.n_id) # p-value = 0.6344


#### REGRESSÃO SIMPLES ####

produt.gini <- lm(c_sid~gini,data = dadosreg_s)
summary(produt.gini) # p-value: 0.5315
residuos_gini <- residuals(produt.gini)
shapiro.test(residuos_gini) #NORMAL

produt.sr <- lm(c_sid~SR,data = dadosreg_s)
summary(produt.sr) # p-value: 0.619
residuos_sr <- residuals(produt.sr)
shapiro.test(residuos_sr) #NORMAL

produt.pse <- lm(c_sid~PSE,data = dadosreg_s)
summary(produt.pse) # p-value: 0.201
residuos_pse <- residuals(produt.pse)
shapiro.test(residuos_pse) #NORMAL

produt.pd <- lm(c_sid~PD,data = dadosreg_s)
summary(produt.pd) # p-value: 0.428
residuos_pd <- residuals(produt.pd)
shapiro.test(residuos_pd) #NORMAL

produt.mpd <- lm(c_sid~MPD,data = dadosreg_s)
summary(produt.mpd) # p-value: 0.178
residuos_mpd <- residuals(produt.mpd)
shapiro.test(residuos_mpd) #NORMAL

produt.sesmpd <- lm(c_sid~SESMPD,data = dadosreg_s)
summary(produt.sesmpd) # p-value: 0.31
residuos_sesmpd <- residuals(produt.sesmpd)
shapiro.test(residuos_sesmpd) #NORMAL

produt.sesmntd <- lm(c_sid~SESMNTD,data = dadosreg_s)
summary(produt.sesmntd) # p-value: 0.717
residuos_sesmntd <- residuals(produt.sesmntd)
shapiro.test(residuos_sesmntd) #NORMAL

produt.ppt <- lm(c_sid~ppt,data = dadosreg_s)
summary(produt.ppt) # p-value: 0.279
residuos_ppt <- residuals(produt.ppt)
shapiro.test(residuos_ppt) #NORMAL
