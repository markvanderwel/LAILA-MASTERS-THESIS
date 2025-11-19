## WOOD DENSITY ###

install.packages("BIOMASS")
library("BIOMASS")

# Input dos dados

teste <- read.csv("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Dados/getwd.csv", sep = ";")

# Compute the Wood Density up to the genus level and give the mean wood density of the dataset

WD2 <- getWoodDensity(
  genus = teste$genus,
  species = teste$species,
  region = "SouthAmericaTrop")

setwd("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Dados")
write.csv(WD, file = "wd_biomass")

mean(x = WD$sdWD) #0.1069175
mean(x = WD2$sdWD) #0.1097415

?getWoodDensity
?BIOMASS

dados <- read.csv("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Dados/TRABALHO EXCEL/comunidade.csv", sep = ";")

# Calcular Biomassa (Mg ou ton)

AGB <- computeAGB(D, WD, H)

agb_chave <- computeAGB(dados$dap,dados$getwd,dados$altura)
agb_chave <- as.matrix(agb_chave)

library(tidyverse)

dados <- dados%>%
  mutate(biomassa_g_c=agb_chave*1000)

write.csv(dados, file = "dados.csv")
