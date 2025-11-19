##############################################
#### DIVERSIDADE FUNCIONAL MESTRADO LAÍLA ####
##############################################

# Carregar pacotes
library(readxl)
library(FD)

setwd("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/Resultados funcional")

#### FUNCIONAL MICO-LEÃO + Exportação para Excel ####

# === 1. Carregar pacotes ===
library(FD)
library(openxlsx)

# === 2. Leitura dos dados ===
comunidade_ml <- read.csv("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/comunidade_ml.csv", 
                          row.names = 1, header = TRUE, sep = ";", check.names = FALSE)

funcional_ml <- read.csv("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/funcional_ml.csv", 
                         row.names = 1, sep = ";")

# === 3. Atributos usados ===
atributos <- c("WD", "SLA", "LDMC")

# === 4. Rodar dbFD combinado ===
resultado_fd_ml <- dbFD(
  x = funcional_ml[, atributos],
  a = comunidade_ml,
  calc.FRic = TRUE,
  calc.FDiv = TRUE,
  calc.CWM  = TRUE,
  m = "max",
  stand.x = TRUE,
  corr = "cailliez"
)

# === 5. Criar planilha Excel e adicionar abas ===
wb <- createWorkbook()

# Adicionar métricas combinadas
addWorksheet(wb, "FDis_Combinado")
writeData(wb, "FDis_Combinado", resultado_fd_ml$FDis)

addWorksheet(wb, "FRic_Combinado")
writeData(wb, "FRic_Combinado", resultado_fd_ml$FRic)

addWorksheet(wb, "FDiv_Combinado")
writeData(wb, "FDiv_Combinado", resultado_fd_ml$FDiv)

addWorksheet(wb, "CWM_por_atributo")
writeData(wb, "CWM_por_atributo", resultado_fd_ml$CWM)

# === 6. Cálculo separado por atributo (sem CWM) ===
for (attr in atributos) {
  cat("Calculando para:", attr, "\n")
  dados_attr <- funcional_ml[, attr, drop = FALSE]
  
  res <- dbFD(x = dados_attr,
              a = comunidade_ml,
              calc.FRic = TRUE,
              calc.FDiv = TRUE,
              calc.CWM = FALSE,
              m = "max",
              stand.x = TRUE,
              corr = "cailliez")
  
  # Adicionar aba com resultados
  addWorksheet(wb, paste0("FDis_", attr))
  writeData(wb, paste0("FDis_", attr), res$FDis)
  
  addWorksheet(wb, paste0("FRic_", attr))
  writeData(wb, paste0("FRic_", attr), res$FRic)
  
  addWorksheet(wb, paste0("FDiv_", attr))
  writeData(wb, paste0("FDiv_", attr), res$FDiv)
}

# === 7. Salvar o arquivo Excel ===
saveWorkbook(wb, file = "Resultados_Funcional_MicoLeao.xlsx", overwrite = TRUE)

###

#### FUNCIONAL REGUA + Exportação para Excel ####

# === 1. Pacotes ===
library(FD)
library(openxlsx)

# === 2. Leitura dos dados ===
comunidade_regua <- read.csv("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/comunidade_regua.csv", 
                             row.names = 1, sep = ";", check.names = FALSE)

funcional_regua <- read.csv("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/funcional_regua.csv", 
                            row.names = 1, sep = ";")

# === 3. Atributos usados ===
atributos <- c("WD", "SLA", "LDMC")

# === 4. Rodar dbFD combinado ===
resultado_fd_combinado <- dbFD(
  x = funcional_regua[, atributos],
  a = comunidade_regua,
  calc.FRic = TRUE,
  calc.FDiv = TRUE,
  calc.CWM  = TRUE,
  m = "max",
  stand.x = TRUE,
  corr = "cailliez"
)

# === 5. Criar planilha Excel e adicionar abas ===
wb <- createWorkbook()

# Adicionar métricas combinadas
addWorksheet(wb, "FDis_Combinado")
writeData(wb, "FDis_Combinado", resultado_fd_combinado$FDis)

addWorksheet(wb, "FRic_Combinado")
writeData(wb, "FRic_Combinado", resultado_fd_combinado$FRic)

addWorksheet(wb, "FDiv_Combinado")
writeData(wb, "FDiv_Combinado", resultado_fd_combinado$FDiv)

addWorksheet(wb, "CWM_por_atributo")
writeData(wb, "CWM_por_atributo", resultado_fd_combinado$CWM)

# === 6. Cálculo separado para cada atributo (sem CWM) ===
for (attr in atributos) {
  cat("Calculando para:", attr, "\n")
  dados_attr <- funcional_regua[, attr, drop = FALSE]
  
  res <- dbFD(x = dados_attr,
              a = comunidade_regua,
              calc.FRic = TRUE,
              calc.FDiv = TRUE,
              calc.CWM = FALSE,
              m = "max",
              stand.x = TRUE,
              corr = "cailliez")
  
  # Adicionar cada resultado como nova aba
  addWorksheet(wb, paste0("FDis_", attr))
  writeData(wb, paste0("FDis_", attr), res$FDis)
  
  addWorksheet(wb, paste0("FRic_", attr))
  writeData(wb, paste0("FRic_", attr), res$FRic)
  
  addWorksheet(wb, paste0("FDiv_", attr))
  writeData(wb, paste0("FDiv_", attr), res$FDiv)
}

# === 7. Salvar arquivo final ===
saveWorkbook(wb, file = "Resultados_Funcional_REGUA.xlsx", overwrite = TRUE)



#########################
## LINEAR MIXED MODELS ##
#########################

library(writexl)
library(readxl)
dadosreg <- read_excel("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/reg_funcional.xlsx")

#### PRODUTIVIDADE ####

sespd.fdis <- lm(SESPD~fdis,data = dadosreg)
summary(sespd.fdis) # p-value: 0.6745  
sespd.fric <- lm(SESPD~fric,data = dadosreg)
summary(sespd.fric) # p-value: 0.7363
sespd.fdiv <- lm(SESPD~fdiv,data = dadosreg)
summary(sespd.fdiv) # p-value: 0.04788 *

produt.fdis <- lm(log_produt~fdis,data = dadosreg)
summary(produt.fdis) # p-value: 0.323     
residuos_fdis <- residuals(produt.fdis)
shapiro.test(residuos_fdis) #NORMAL

produt.fric <- lm(log_produt~fric,data = dadosreg)
summary(produt.fric) # p-value: 0.878 
residuos_fric <- residuals(produt.fric)
shapiro.test(residuos_fric) #NORMAL

produt.fdiv <- lm(log_produt~fdiv,data = dadosreg)
summary(produt.fdiv) # p-value: 0.509
residuos_fdiv <- residuals(produt.fdiv)
shapiro.test(residuos_fdiv) #NORMAL

produt.cwm_wd <- lm(log_produt~cwm_wd,data = dadosreg)
summary(produt.cwm_wd) # p-value: 0.167 
residuos_cwm_wd <- residuals(produt.cwm_wd)
shapiro.test(residuos_cwm_wd) #NORMAL

produt.cwm_sla <- lm(log_produt~cwm_sla,data = dadosreg)
summary(produt.cwm_sla) # p-value: 0.579 
residuos_cwm_sla <- residuals(produt.cwm_sla)
shapiro.test(residuos_cwm_sla) #NORMAL

produt.cwm_ldmc <- lm(log_produt~cwm_ldmc,data = dadosreg)
summary(produt.cwm_ldmc) # p-value: 0.146
residuos_cwm_ldmc <- residuals(produt.cwm_ldmc)
shapiro.test(residuos_cwm_ldmc) #NORMAL

produt.fdis_wd <- lm(log_produt~fdis_wd,data = dadosreg)
summary(produt.fdis_wd) # p-value: 0.699
residuos_fdis_wd <- residuals(produt.fdis_wd)
shapiro.test(residuos_fdis_wd) #NORMAL

produt.fric_wd <- lm(log_produt~fric_wd,data = dadosreg)
summary(produt.fric_wd) # p-value: 0.806
residuos_fric_wd <- residuals(produt.fric_wd)
shapiro.test(residuos_fric_wd) #NORMAL

produt.fdis_sla <- lm(log_produt~fdis_sla,data = dadosreg)
summary(produt.fdis_sla) # p-value: 0.214
residuos_fdis_sla <- residuals(produt.fdis_sla)
shapiro.test(residuos_fdis_sla) #NORMAL

produt.fric_sla <- lm(log_produt~fric_sla,data = dadosreg)
summary(produt.fric_sla) # p-value: 0.333 
residuos_fric_sla <- residuals(produt.fric_sla)
shapiro.test(residuos_fric_sla) #NORMAL

produt.fdis_ldmc <- lm(log_produt~fdis_ldmc,data = dadosreg)
summary(produt.fdis_ldmc) # p-value: 0.0483 *
residuos_fdis_ldmc <- residuals(produt.fdis_ldmc)
shapiro.test(residuos_fdis_ldmc) #NORMAL

produt.fric_ldmc <- lm(log_produt~fric_ldmc,data = dadosreg)
summary(produt.fric_ldmc) # p-value: 0.0301 * 
residuos_fric_ldmc <- residuals(produt.fric_ldmc)
shapiro.test(residuos_fric_ldmc) #NORMAL

produt.SR <- lm(log_produt~SR,data = dadosreg)
summary(produt.SR) # p-value: 0.7908
residuos_sr <- residuals(produt.SR)
shapiro.test(residuos_sr) #NORMAL

t <- lm(log_produt~cwm_ldmc+cwm_wd,data = dadosreg)
summary(t) # p-value: 0.044

t2 <- lm(log_produt ~ cwm_ldmc + cwm_wd + fdis_ldmc + fric_ldmc, data = dadosreg)
anova(t, t2)        # Likelihood-ratio test
AIC(t); AIC(t2)     # Compare AIC (should drop if FD adds value)
performance::r2(t); performance::r2(t2)   # ΔR²

# 1.1 FDis only
m_fdis <- lm(log_produt ~ cwm_ldmc + cwm_wd + fdis_ldmc, data = dadosreg)
summary(m_fdis)
AIC(m_fdis)
performance::r2(m_fdis)

# 1.2 FRic only
m_fric <- lm(log_produt ~ cwm_ldmc + cwm_wd + fric_ldmc, data = dadosreg)
summary(m_fric)
AIC(m_fric)
performance::r2(m_fric)

# 2. Nonlinear (quadratic) test for FDis
m_fdis_quad <- lm(log_produt ~ cwm_ldmc + cwm_wd +
                    fdis_ldmc + I(fdis_ldmc^2), data = dadosreg)
summary(m_fdis_quad)
anova(t, m_fdis_quad)  # compare to CWM-only baseline


##

t  <- lm(log_produt ~ cwm_ldmc + cwm_wd, data = dadosreg)
t2 <- lm(log_produt ~ fdis_ldmc, data = dadosreg)
AIC(t, t2)
summary(t)


# site as categorical variable

m_site <- lm(log_produt ~ cwm_ldmc + cwm_wd + sitemis, data = dadosreg)
summary(m_site)


m_fd_site <- lm(log_produt ~ cwm_ldmc + cwm_wd + fdis_ldmc * sitemis, data = dadosreg)
anova(m_site, m_fd_site)

by(dadosreg, dadosreg$sitemis, function(df) summary(lm(log_produt ~ cwm_ldmc + cwm_wd + fdis_ldmc, data = df)))

### Plot ###

library(ggplot2)

ggplot(dadosreg, aes(x = cwm_ldmc, y = log_produt, color = sitemis)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_classic() +
  labs(x = "CWM LDMC", y = "Log(Biomass/age)", color = "Site")

ggplot(dadosreg, aes(x = fdis_ldmc, y = log_produt, color = sitemis)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_classic() +
  labs(x = "FDis LDMC", y = "Log(Biomass/age)", color = "Site")


library(ggplot2)
library(patchwork)   # to combine plots side by side

# Panel A – Mass-ratio (CWM LDMC)
p1 <- ggplot(dadosreg, aes(x = cwm_ldmc, y = log_produt, color = sitemis)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "CWM LDMC", y = "Log(Biomass/age)", color = "Site",
       title = "(a) Mass-ratio effect (CWM)") +
  theme_classic() +
  theme(legend.position = "top")

# Panel B – Complementarity (FDis LDMC)
p2 <- ggplot(dadosreg, aes(x = fdis_ldmc, y = log_produt, color = sitemis)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "FDis LDMC", y = NULL, color = "Site",
       title = "(b) Complementarity effect (FD)") +
  theme_classic() +
  theme(legend.position = "none")   # shared legend from p1

# Combine
p1 + p2

####### Effect plots for each variable #######

library(ggplot2)
library(ggeffects)

# Get partial predictions (controlling for the other variable)
eff_ldmc <- ggpredict(t, terms = "cwm_ldmc")
eff_wd   <- ggpredict(t, terms = "cwm_wd")

# Plot both effects side by side
p1 <- ggplot(eff_ldmc, aes(x, predicted)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  labs(x = "CWM LDMC", y = "Predicted log(Biomass/age)",
       title = "(a) Effect of CWM LDMC") +
  theme_classic()

p2 <- ggplot(eff_wd, aes(x, predicted)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  labs(x = "CWM WD", y = "Predicted log(Biomass/age)",
       title = "(b) Effect of CWM WD") +
  theme_classic()

library(patchwork)
p1 + p2


#######

install.packages("plotly")
library(plotly)
library(tidyverse)

library(plotly)

# 1) Fit (you already have it)
t <- lm(log_produt ~ cwm_ldmc + cwm_wd, data = dadosreg)

library(plotly)


# prediction grid
xseq <- seq(min(dadosreg$cwm_ldmc, na.rm=TRUE),
            max(dadosreg$cwm_ldmc, na.rm=TRUE), length.out=40)
yseq <- seq(min(dadosreg$cwm_wd,   na.rm=TRUE),
            max(dadosreg$cwm_wd,   na.rm=TRUE), length.out=40)

grid  <- expand.grid(cwm_ldmc = xseq, cwm_wd = yseq)
zpred <- matrix(predict(t, newdata = grid),
                nrow = length(xseq), ncol = length(yseq))

# plot
p <- plot_ly() |>
  add_markers(data = dadosreg,
              x = ~cwm_ldmc, y = ~cwm_wd, z = ~log_produt,
              color = ~sitemis, marker = list(size = 4, opacity = 0.9),
              name = "Observed") |>
  add_trace(x = xseq, y = yseq, z = zpred,
            type = "surface", opacity = 0.6, showscale = FALSE,
            name = "LM surface") |>
  plotly::layout(scene = list(
    xaxis = list(title = "CWM LDMC"),
    yaxis = list(title = "CWM WD"),
    zaxis = list(title = "log(Biomass/age)")
  ))

install.packages("webshot2")
library(webshot2)

setwd("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/Resultados funcional/Gráficos")
htmlwidgets::saveWidget(p, "3D_model.html")   # salva o gráfico interativo
webshot("3D_model.html", "3D_model.jpeg", vwidth = 1600, vheight = 1200)



### FACET WRAP ###
library(tidyverse)
library(ggplot2)
setwd("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/Resultados funcional/Gráficos")

# SR
dadosreg %>%
  ggplot(aes(SR,log_produt))+
  geom_point(size=3,alpha=0.5)+ #alpha é transparência
  geom_smooth(method = lm,se=F)+ #linear model, se = F não tem intervalo de confi
  facet_wrap(~sitemis)+ # divide em dois gráficos diferentes
  labs(x="SR",
       y="Productivity",
       title = "Species Richness and Productivity in ecological restoration plantings")+
  theme_bw()

ggsave("SR_prod_facetwrsp.jpg", width = 15, height = 10, units = "cm", dpi = 300)

# FDIS_SLA
dadosreg %>%
  ggplot(aes(fric_wd,log_produt))+
  geom_point(size=3,alpha=0.5)+ #alpha é transparência
  geom_smooth(method = lm,se=F)+ #linear model, se = F não tem intervalo de confi
  facet_wrap(~sitemis)+ # divide em dois gráficos diferentes
  labs(x="SR",
       y="Productivity",
       title = "Functional Richness - SLA")+
  theme_bw()

ggsave("SR_prod_facetwrsp.jpg", width = 15, height = 10, units = "cm", dpi = 300)


# Separar por sítio e rodar regressões
modelo_ml <- lm(log_produt ~ fric_wd, data = dadosreg %>% filter(sitemis == "ML"))
modelo_regua <- lm(log_produt ~ fric_wd, data = dadosreg %>% filter(sitemis == "REGUA"))

# Ver os resumos
summary(modelo_ml)
summary(modelo_regua)

# FRic_LDMC
dadosreg %>%
  ggplot(aes(fric_ldmc,log_produt))+
  geom_point(size=3,alpha=0.5)+ #alpha é transparência
  geom_smooth(method = lm,se=F)+ #linear model, se = F não tem intervalo de confi
  facet_wrap(~sitemis)+ # divide em dois gráficos diferentes
  labs(x="Fric",
       y="Productivity",
       title = "Functional Richness and Productivity in ecological restoration plantings")+
  theme_bw()

ggsave("FRIC_LDMC_facetwrap.jpg", width = 15, height = 10, units = "cm", dpi = 300)

# Separar por sítio e rodar regressões
modelo_ml <- lm(log_produt ~ fric_ldmc, data = dadosreg %>% filter(sitemis == "ML"))
modelo_regua <- lm(log_produt ~ fric_ldmc, data = dadosreg %>% filter(sitemis == "REGUA"))

# Ver os resumos
summary(modelo_ml)
summary(modelo_regua)

### GRÁFICO NORMAL

# FDis
dadosreg %>%
  ggplot(aes(fdis_ldmc,log_produt))+
  geom_point(size=3,alpha=0.5)+ 
  geom_smooth(method = lm,se=T,colour="#7FFFD4")+
  labs(x="FDis - Leaf Dry Matter Content (LDMC)",
       y="log Productivity (g/m²/ano)",
       title = "The Relationship between Functional Diversity and Productivity")+
  annotate("text",x=1,y=4.58,label="R² = 0.17 
p-valor < 0.05")+
  theme_replace()

ggsave("FDis_LDMC.png",width = 15,height = 10,units = "cm")


# FRic LDMC
dadosreg %>%
  ggplot(aes(fric_ldmc,log_produt))+
  geom_point(size=3,alpha=0.5)+ 
  geom_smooth(method = lm,se=T,colour="yellow")+
  labs(x="FRic - Leaf Dry Matter Content (LDMC)",
       y="log Productivity (g/m²/ano)",
       title = "The Relationship between Functional Richness and Productivity")+
  annotate("text",x=2.75,y=4.58,label="R² = 0.20 
p-valor = 0.03")+
  theme_replace()

ggsave("FRic_LDMC.png",width = 15,height = 10,units = "cm")

# CWM LDMC
dadosreg %>%
  ggplot(aes(cwm_ldmc,log_produt))+
  geom_point(size=3,alpha=0.5)+ 
  geom_smooth(method = lm,se=T,colour="red")+
  labs(x="FDis - Leaf Dry Matter Content (LDMC)",
       y="log Productivity (g/m²/ano)",
       title = "The Relationship between Functional Diversity and Productivity")+
  annotate("text",x=1,y=4.58,label="R² = 0.17 
p-valor < 0.05")+
  theme_replace()

###########################

setwd("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/Resultados funcional/Gráficos")

# Ajuste do modelo
fit <- lm(log_produt ~ fric_ldmc, data = dadosreg)
summary(fit)

# Extraindo R² e p-valor
r2 <- summary(fit)$r.squared
p  <- coef(summary(fit))["fric_ldmc", "Pr(>|t|)"]

# Texto formatado
subtxt <- sprintf("R² = %.2f\np = %s",
                  r2,
                  ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

# plot
g <- dadosreg %>%
  ggplot(aes(fric_ldmc, log_produt)) +
  geom_point(size = 3, alpha = 0.5) + 
  geom_smooth(method = "lm", se = TRUE, colour = "yellow") +
  labs(
    x = "FRic - Leaf Dry Matter Content (LDMC)",
    y = expression("log Net Biomass (g m"^-2*" year"^-1*")"),
    title = "The Relationship between Functional Richness and Productivity"
  ) +
  annotate("text", 
           x = max(dadosreg$fric_ldmc, na.rm = TRUE) * 0.45, 
           y = max(dadosreg$log_produt, na.rm = TRUE) * 1.002, 
           label = subtxt, 
           hjust = 1, vjust = 1,
           size = 3.5) +
  theme_minimal(base_size = 11)

g

# Salvar
ggsave("FRic_LDMC.jpeg", g, width = 15, height = 10, units = "cm", dpi = 600)


### FDIS ###

fit <- lm(log_produt ~ fdis_ldmc, data = dadosreg)
summary(fit)

# Extraindo R² e p-valor
r2 <- summary(fit)$r.squared
p  <- coef(summary(fit))["fdis_ldmc", "Pr(>|t|)"]

# Texto formatado
subtxt <- sprintf("R² = %.2f\np = %s",
                  r2,
                  ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

# Gráfico
g <- dadosreg %>%
  ggplot(aes(fdis_ldmc, log_produt)) +
  geom_point(size = 3, alpha = 0.5) + 
  geom_smooth(method = "lm", se = TRUE, colour = "brown") +
  labs(
    x = "FDis - Leaf Dry Matter Content (LDMC)",
    y = expression("log Net Biomass (g m"^-2*" year"^-1*")"),
    title = "The Relationship between Functional Richness and Productivity"
  ) +
  annotate("text", 
           x = max(dadosreg$fric_ldmc, na.rm = TRUE) * 0.1, 
           y = max(dadosreg$log_produt, na.rm = TRUE) * 1.002, 
           label = subtxt, 
           hjust = 1, vjust = 1,
           size = 3.5) +
  theme_minimal(base_size = 11)

g

# Salvar
ggsave("FDis_LDMC.jpeg", g, width = 15, height = 10, units = "cm", dpi = 600)


## CWM LDMC ##

# 1) Fit the model
fit <- lm(log_produt ~ cwm_ldmc, data = dadosreg)
sm  <- summary(fit)

# 2) Extract stats safely
r2    <- sm$r.squared
adjr2 <- sm$adj.r.squared
p     <- sm$coefficients["cwm_ldmc", "Pr(>|t|)"]

subtxt <- sprintf("R² = %.2f | adj. R² = %.2f\np = %s",
                  r2, adjr2,
                  ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

# 3) Nice positions for the annotation inside the panel
x_annot <- min(dadosreg$cwm_ldmc,  na.rm = TRUE)
y_annot <- max(dadosreg$log_produt, na.rm = TRUE)

# 4) Plot
g <- dadosreg %>%
  ggplot(aes(x = cwm_ldmc, y = log_produt)) +
  geom_point(size = 3, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.9, colour = "goldenrod") +
  labs(
    x = "CWM — Leaf Dry Matter Content (LDMC)",
    y = expression("log Net Biomass (g m"^-2*" year"^-1*")"),
    title = "Relationship between CWM-LDMC and Productivity"
  ) +
  annotate("label",
           x = x_annot, y = y_annot,
           hjust = 0, vjust = 1,
           label = subtxt, size = 3.5) +
  theme_minimal(base_size = 12)

g



######################################
### PHYLOGENETIC SIGNAL FOR TRAITS ###
######################################

## MONTAR A ÁRVORE SÓ COM AS ESPÉCIES DE ML E REGUA

# load the packages

library("V.PhyloMaker2")
library(writexl)

# input the sample species list
example <- read.csv("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/filogenia_ml_regua.csv")

### generate a phylogeny for the sample species list
tree <- phylo.maker(example, tree = GBOTB.extended.TPL,output.sp.list = TRUE,nodes = nodes.info.1.TPL, scenarios="S3")
tree_ok <- tree$scenario.3
summary(tree_ok)

setwd("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística")
dados <- read.csv("funcional_ml_regua.csv", row.names = 1,header = T, sep = ";")

library(ape)

tree_teste <- multi2di(tree_ok)

library(picante)

resultado <- multiPhylosignal(dados[, c("WD", "SLA", "LDMC")], tree_teste)
# LDMC - SINAL FORTE

#
#
#
### DECOUPLED EFFECTS ###
#
#
#

# Load required packages
library(devtools)
devtools::install_github("carmonalab/decouple")

library(decouple)
library(ape)        # For phylogenetic tree manipulation
library(vegan)      # For dissimilarity calculations
library(picante)    # For MPD calculation
library(lme4)       # For mixed models
library(MuMIn)      # For AICc and model R²

# Prepare input data

# Extract only LDMC trait from the trait matrix
traits_ldmc <- as.matrix(traits_df["LDMC"])  # species x 1 trait

# Ensure community matrix is species x plots
abund_matrix <- comunidade_ldmc

# Ensure species names match across tree, traits, and community matrix
all(rownames(traits_ldmc) %in% tree_teste$tip.label)  # Should return TRUE
all(rownames(abund_matrix) %in% tree_teste$tip.label) # Should also return TRUE

# Run the decoupling function
res <- decouple(
  traits = traits_ldmc,
  phylo = tree_teste,
  abund = abund_matrix,
  dist_method = "euclidean",  # for single continuous trait
  stand.x = TRUE,             # standardize trait
  stand.P = TRUE              # standardize phylogenetic distances
)

# Extract decoupled dissimilarity matrices
dcFdist <- res$dcFdist  # Functional dissimilarity independent of phylogeny
dcPdist <- res$dcPdist  # Phylogenetic dissimilarity independent of traits

# Calculate mean pairwise dissimilarity (MPD) per plot
# transpose community matrix to plots x species as required by picante
mpd_dcF <- mpd(samp = t(abund_matrix), dis = dcFdist)
mpd_dcP <- mpd(samp = t(abund_matrix), dis = dcPdist)

# Add results to your data frame for modeling
dadosreg$mpd_dcF <- mpd_dcF[rownames(dadosreg)]
dadosreg$mpd_dcP <- mpd_dcP[rownames(dadosreg)]

# Fit mixed models using decoupled metrics
modelo_dcF <- lmer(log_produt ~ mpd_dcF + (1 | sitemis), data = dadosreg, REML = FALSE)
modelo_dcP <- lmer(log_produt ~ mpd_dcP + (1 | sitemis), data = dadosreg, REML = FALSE)

# Summarize results
summary(modelo_dcF)
summary(modelo_dcP)

# Compare models
r.squaredGLMM(modelo_dcF)
r.squaredGLMM(modelo_dcP)
AICc(modelo_dcF)
AICc(modelo_dcP)



################################################################################


### COEFPLOT ###

library(dplyr)
library(ggplot2)
library(purrr)
library(stringr)

# ---------- 1) LISTA DE PREDITORES ----------
preds <- c("SR",
           "fdis","fric","fdiv",
           "cwm_wd","cwm_sla","cwm_ldmc",
           "fdis_wd","fric_wd",
           "fdis_sla","fric_sla",
           "fdis_ldmc","fric_ldmc")

# Rótulos amigáveis para o gráfico
pretty_lab <- c(
  SR           = "Species richness (SR)",
  fdis         = "FDis (all traits)",
  fric         = "FRic (all traits)",
  fdiv         = "FDiv (all traits)",
  cwm_wd       = "CWM — Wood density",
  cwm_sla      = "CWM — SLA",
  cwm_ldmc     = "CWM — LDMC",
  fdis_wd      = "FDis — WD",
  fric_wd      = "FRic — WD",
  fdis_sla     = "FDis — SLA",
  fric_sla     = "FRic — SLA",
  fdis_ldmc    = "FDis — LDMC",
  fric_ldmc    = "FRic — LDMC"
)

# ---------- 2) FUNÇÃO PARA AJUSTAR E EXTRAIR COEF ----------
fit_and_tidy <- function(var, dat) {
  form <- as.formula(paste0("log_produt ~ ", var))
  fit  <- lm(form, data = dat)
  sm   <- summary(fit)
  ci   <- suppressWarnings(confint(fit, level = 0.95))
  # Pega só o termo do preditor
  est  <- sm$coefficients[var, "Estimate"]
  se   <- sm$coefficients[var, "Std. Error"]
  pval <- sm$coefficients[var, "Pr(>|t|)"]
  lo   <- ci[var, 1]
  hi   <- ci[var, 2]
  tibble(
    predictor = var,
    estimate  = est,
    conf.low  = lo,
    conf.high = hi,
    p.value   = pval
  )
}

# ---------- 3) TABELA DE COEFICIENTES ----------
coef_tab <- map_dfr(preds, fit_and_tidy, dat = dadosreg) %>%
  mutate(label      = pretty_lab[predictor],
         is_ldmc    = str_detect(predictor, "ldmc"),
         sig        = p.value < 0.05,
         label_plot = ifelse(sig, paste0("**", label, "**"), label))

# Ordena por magnitude do efeito (absoluta) ou como preferir
coef_tab <- coef_tab %>%
  arrange(estimate) %>%
  mutate(label_plot = factor(label_plot, levels = label_plot))

# ---------- 4) PLOT ----------
g_coef <- ggplot(coef_tab,
                 aes(x = estimate, y = label_plot)) +
  geom_vline(xintercept = 0, linetype = 2, color = "gray70") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.2,
                 alpha  = ifelse(coef_tab$sig, 1, 0.5),
                 color  = ifelse(coef_tab$is_ldmc, "#E5C100", "gray55")) +
  geom_point(size = 2.8,
             alpha = ifelse(coef_tab$sig, 1, 0.7),
             color = ifelse(coef_tab$is_ldmc, "#E5C100", "gray25")) +
  labs(x = "Slope (β) ± 95% CI",
       y = NULL,
       title = "Effect sizes of trait metrics on productivity",
       subtitle = "Highlighted in yellow: LDMC metrics; Bold labels: p < 0.05") +
  theme_minimal(base_size = 12) +
  theme(
    axis.title.x   = element_text(face = "bold"),
    plot.subtitle  = element_text(color = "gray30"),
    panel.grid.minor = element_blank()
  )

g_coef
ggsave("coefplot_produt_raw.jpeg", g_coef, width = 16, height = 11, units = "cm", dpi = 600)




####### Not used anymore ############




### MODELO MISTO # SITE AS RANDOM VARIABLE ###

library(lmerTest) #MODELO MISTO
library(MuMIn) #AICc

# Modelo nulo sem os preditores fixos
modelo_nulo <- lmer(log_produt ~ (1 | sitemis), data = dadosreg, REML = FALSE)
summary(modelo_nulo)
AICc(modelo_nulo) # -12.20405


misto1 <- lmer(log_produt ~ fric_wd + (1 | sitemis), data = dadosreg, REML = FALSE)
summary(misto1)  # p-value = 0.68
r.squaredGLMM(misto1) # 0.5769652
AICc(misto1) # -9.418423

misto2 <- lmer(log_produt ~ fric_sla + (1 | sitemis), data = dadosreg, REML = FALSE)
summary(misto2)  # p-value = 0.816
r.squaredGLMM(misto2) # 0.5712046
AICc(misto2) # -9.296783

misto3 <- lmer(log_produt ~ fric_ldmc + (1 | sitemis), data = dadosreg, REML = FALSE)
summary(misto3)  # p-value = 0.0643 .
r.squaredGLMM(misto3) # 0.5857842
AICc(misto3) # -12.76521

misto4 <- lmer(log_produt ~ fdis_wd + (1 | sitemis), data = dadosreg, REML = FALSE)
summary(misto4)  # p-value = 0.839    
r.squaredGLMM(misto4) # 0.5660719
AICc(misto4) # -9.287176

misto5 <- lmer(log_produt ~ fdis_sla + (1 | sitemis), data = dadosreg, REML = FALSE)
summary(misto5)  # p-value = 0.617        
r.squaredGLMM(misto5) # 0.5658284
AICc(misto5) # -9.291259

misto6 <- lmer(log_produt ~ fdis_ldmc + (1 | sitemis), data = dadosreg, REML = FALSE)
summary(misto6)  # p-value = 0.229            
r.squaredGLMM(misto6) # 0.5517717
AICc(misto6) # -10.71563

misto7 <- lmer(log_produt ~ cwm_wd + (1 | sitemis), data = dadosreg, REML = FALSE)
summary(misto7)  # p-value = 0.658                
r.squaredGLMM(misto7) # 0.5550451
AICc(misto7) # -9.443239

misto8 <- lmer(log_produt ~ cwm_sla + (1 | sitemis), data = dadosreg, REML = FALSE)
summary(misto8)  # p-value = 0.592                
r.squaredGLMM(misto8) # 0.5101121
AICc(misto8) # -10.33825

misto9 <- lmer(log_produt ~ cwm_ldmc + (1 | sitemis), data = dadosreg, REML = FALSE)
summary(misto9)  # p-value = 0.955                    
r.squaredGLMM(misto9) # 0.5641951
AICc(misto9) # -9.248183

misto10 <- lmer(log_produt ~ fdis + (1 | sitemis), data = dadosreg, REML = FALSE)
summary(misto10)  # p-value = 0.635                        
r.squaredGLMM(misto10) # 0.5665168
AICc(misto10) # -9.32659

misto11 <- lmer(log_produt ~ fric + (1 | sitemis), data = dadosreg, REML = FALSE)
summary(misto11)  # p-value = 0.280460                            
r.squaredGLMM(misto11) # 0.606324
AICc(misto11) # -10.66883

misto12 <- lmer(log_produt ~ fdiv + (1 | sitemis), data = dadosreg, REML = FALSE)
summary(misto12)  # p-value = 0.991                                
r.squaredGLMM(misto12) # 0.562682
AICc(misto12) # -9.353998

misto13 <- lmer(log_produt ~ cwm_sla + fric_ldmc + (1 | sitemis), data = dadosreg, REML = FALSE)
summary(misto13)
anova(modelo_nulo, misto13) # p-value = 0.1546
r.squaredGLMM(misto13) # 0.5101121
AICc(misto13) # -10.33825

misto14 <- lmer(log_produt ~ SR + (1 | sitemis), data = dadosreg, REML = FALSE)
summary(misto14) # p-value = 0.143
r.squaredGLMM(misto14) # 0.6546653
AICc(misto14) # -11.37825

misto15 <- lmer(log_produt ~ SESPD + SR + fric_ldmc + (1 | sitemis), data = dadosreg, REML = FALSE)
summary(misto15)
anova(modelo_nulo, misto15) # p-value = 0.1354
r.squaredGLMM(misto15) # 0.6546653
AICc(misto15) # -9.937358

misto15 <- lmer(log_produt ~ PSEab + fric_ldmc + (1 | sitemis), data = dadosreg, REML = FALSE)
summary(misto15)
anova(modelo_nulo, misto15) # p-value = 0.1354
r.squaredGLMM(misto15) # 0.6546653
AICc(misto15) # -9.937358

misto15 <- lmer(log_produt ~ dcP + (1 | sitemis), data = dadosreg, REML = FALSE)
summary(misto15)
anova(modelo_nulo, misto15) # p-value = 0.1779
r.squaredGLMM(misto15) # 0.6529036
AICc(misto15) # -11.0603


