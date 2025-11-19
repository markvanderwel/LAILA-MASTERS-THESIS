## ORGANIZAR ##

# OFICIAL ### PRODUTIVIDADE ### site as random variable

# Modelo nulo sem os preditores fixos
modelo_nulo <- lmer(log_produt ~ (1 | site), data = dadosmisto, REML = FALSE)
summary(modelo_nulo)
AICc(modelo_nulo)

misto1 <- lmer(log_produt ~ silte + valor_t + SR + (1 | site), data = dadosmisto, REML = FALSE)
summary(misto1) # só SR significativo
anova(modelo_nulo, misto) # p-value = 8.162e-07 ***
r.squaredGLMM(misto) # 0.5871057
AICc(misto) # 67.02762

misto2 <- lmer(log_produt ~ silte + SR + (1 | site), data = dadosmisto, REML = FALSE)
summary(misto2)  
anova(modelo_nulo, misto2) # p-value = 0.002972 **
r.squaredGLMM(misto2) # 0.423191
AICc(misto2) # 84.09125


misto3 <- lmer(log_produt ~ SR + PC1_clima + (1 | site), data = dadosmisto, REML = FALSE)
summary(misto3)  # clima não é significativo
anova(modelo_nulo, misto3) # p-value = 0.01057 *
r.squaredGLMM(misto3) # 0.4758896
AICc(misto3) # 86.62853

misto4 <- lmer(log_produt ~ SR + pcps2 + (1 | site), data = dadosmisto, REML = FALSE)
summary(misto4) 
anova(modelo_nulo, misto4) # p-value = 0.001196 **
r.squaredGLMM(misto4) # 0.6224514
AICc(misto4) # 82.27149

misto5 <- lmer(log_produt ~ SR + pcps2 + silte + (1 | site), data = dadosmisto, REML = FALSE)
summary(misto5)
anova(modelo_nulo, misto5) # p-value = 0.0007439 ***
r.squaredGLMM(misto5) # 0.5142652
AICc(misto5) # 81.21954

misto6 <- lmer(log_produt ~ shannon + pcps2 + silte + (1 | site), data = dadosmisto, REML = FALSE)
summary(misto6)
anova(modelo_nulo, misto6) # p-value = 0.008743 **
r.squaredGLMM(misto6) # 0.4556971
AICc(misto6) # 86.47613

#####################

# Carregar pacotes necessários
library(lme4)
library(lmerTest)
library(ggplot2)
library(performance) # Para diagnóstico do modelo

# Ajuste dos modelos
modelo_linear <- lmer(log_produt ~ pcps2 + silte + SR + (1 | site), data = dadosmisto_std)
modelo_quadratico <- lmer(log_produt ~ poly(pcps2, 2) + silte + SR + (1 | site), data = dadosmisto_std)

# Comparação dos modelos
anova(modelo_linear, modelo_quadratico) # Teste de verossimilhança

# Diagnóstico dos resíduos
check_model(modelo_quadratico) # Avaliação gráfica de suposições

# Visualização do modelo linear
ggplot(dadosmisto_std, aes(x = pcps2, y = log_produt)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, color = "red") +
  labs(title = "Relação linear entre PCPS2 e log_produt")

# Visualização do modelo quadrático
ggplot(dadosmisto_std, aes(x = pcps2, y = log_produt)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = "blue") +
  labs(title = "Relação quadrática entre PCPS2 e log_produt")
























# Definindo grupos filogenéticos
grupos <- example$family

# Rodando PCPS
res <- pcps(comunidade, pdis.ord)
summary(res)

# Exportando eixos principais
dadosmisto$pcps1ab <- res$vectors[, 1]
dadosmisto$pcps2ab <- res$vectors[, 2]
write_xlsx(dadosmisto, "dadosmisto.xlsx")

# Calculando variância explicada
eig <- res$values$Eigenvalue
var_exp <- eig / sum(eig)
pcps1_var <- round(var_exp[1] * 100, 2)
pcps2_var <- round(var_exp[2] * 100, 2)

# Scores para sites e espécies
scores_pcps <- scores.pcps(res, choice = c(1,2))
scores_sites <- scores_pcps$scores.sites
scores_species <- scores_pcps$scores.species

# Setas de grupos (famílias)
apg_arrows <- apply(scores_species, 2, function(x) tapply(x, list(grupos), mean)) %>%
  as.data.frame()
apg_arrows$angle <- atan2(apg_arrows$pcps.1, apg_arrows$pcps.2) * (180/pi)
apg_arrows$labels <- c("Anac", "Anno", "Apocy", "Aster", "Bigno", "Bixa", "Borag", "Canna", 
                       "Caric", "Euphorb", "Faba", "Lami", "Laura", "Lecyt", "Malva", 
                       "Melast", "Melia", "Myrta", "Petive", "Polygo", "Primula", "Rubi", 
                       "Ruta", "Sapind", "Solan", "Urtic", "Verbe")

# Juntando os dados em um único data.frame
df_plot <- cbind(as.data.frame(scores_sites),
                 log_produt = dadosmisto$log_produt,
                 SESPDab = dadosmisto$SESPDab)

# Agora sim o ggplot funciona 100%
ggplot(data = df_plot, aes(x = pcps.1, y = pcps.2)) +
  geom_point(aes(color = SESPDab, size = log_produt)) +
  scale_colour_gradient(low = "#0291cb", high = "#ea064f") +
  scale_size_continuous(name = "Biomass (kg)", range = c(2, 8), trans = "log10") +
  labs(x = paste0("PCPS1 (", pcps1_var, "%)"), 
       y = paste0("PCPS2 (", pcps2_var, "%)"),
       title = "Phylogenetic composition of the study areas") +
  geom_text(data = apg_arrows,
            aes(x = pcps.1, y = pcps.2, label = labels), 
            fontface = "bold", color = "black", size = 1) +
  geom_label(data = apg_arrows,
             aes(x = pcps.1, y = pcps.2, label = labels), 
             fontface = "bold", fill = "grey", alpha = 0.5,
             color = "black", size = 3.5) +
  theme_minimal() +
  theme(axis.line = element_line(color = "black", linewidth = 0.8),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 13))










