library(readxl)
library(car)
library(ggeffects)
library(ggplot2)
library(lme4)

dadosmisto <- read_excel("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/dadosmisto.xlsx")

misto.Produt.PSE.ppt.SR.Site <- lmer(log_produt ~ PSE + ppt + SR + (1 | site), data = dadosmisto, REML = FALSE)

# Plotar resíduos do modelo
plot(misto.Produt.PSE.ppt.SR.Site)

# Gerar efeitos marginais para a variável PSE
effects_PSE <- ggpredict(misto.Produt.PSE.ppt.SR.Site, terms = "PSE")
effects_SR <- ggpredict(misto.Produt.PSE.ppt.SR.Site, terms = "SR")
effects_ppt <- ggpredict(misto.Produt.PSE.ppt.SR.Site, terms = "ppt")

# Plotar os efeitos marginais com alterações no título e nos eixos
grafico <- plot(effects_SR) 
setwd("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/graficos")
ggsave("grafico_Predicted.PSE.SR.jpeg", plot = grafico, device = "jpeg", width = 8, height = 6, dpi = 300)

# PSE
grafico <- plot(effects_PSE) +
  ggtitle("Predicted values of productivity (best LMM model)") +
  xlab("PSE") +
  ylab("Productivity (g/m²/year)")

ggsave("Predicted_PSE.jpeg", plot = grafico, device = "jpeg", width = 8, height = 6, dpi = 300)

##### SR
grafico <- plot(effects_SR) +
  ggtitle("Predicted values of productivity (best LMM model)") +
  xlab("SR") +
  ylab("Productivity (g/m²/year)")


ggsave("Predicted_SR.jpeg", plot = grafico, device = "jpeg", width = 8, height = 6, dpi = 300)


##### SR
grafico <- plot(effects_ppt) +
  ggtitle("Predicted values of productivity (best LMM model)") +
  xlab("ppt") +
  ylab("Productivity (g/m²/year)")


ggsave("Predicted_ppt.jpeg", plot = grafico, device = "jpeg", width = 8, height = 6, dpi = 300)


# Criar gráficos de resíduos parciais para as variáveis do modelo

jpeg("crPlots_M6_ppt_dec.jpeg", width = 800, height = 600)

crPlots(M6.ppt.dec)

dev.off()

