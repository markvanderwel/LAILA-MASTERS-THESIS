### Biomassa Chave ###

tree_data <- read.csv("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Dados/Biomassa_Chave.csv", header = T, sep = ";")

library(BIOMASS)

agb <- computeAGB(D = tree_data$DAP, WD = tree_data$WD, H = tree_data$H, coord = NULL)

tree_data$AGB <- agb*1000