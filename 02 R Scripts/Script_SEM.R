###########
### SEM ###
###########

install.packages("tidySEM")
install.packages("semTools")
installed.packages("semPlot")
install.packages("qgraph")

library(readxl)
library(readr)
library(ggplot2)
library(dplyr)

library(lavaan)
library(tidySEM)
library(semTools)
library(semPlot)
library(qgraph)

alldata <- read_excel("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/dadosmisto.xlsx")

# Selecting relevant variables by their name in the dataset itself
variable_list <- c("log_produt", "SR", "ppt", "faba", "SESPD", "gini","valor_s")
# Subsetting the dataframe
df <- alldata[variable_list]
# Renaming columns to make them easily referenced for SEM analysis
colnames(df) <- c("produt", "sr", "ppt", "faba", "pd", "gini","v_s")

shapiro.test(df$sr)
hist(df$sr)
shapiro.test(df$ppt)
hist(df$ppt) # não é normal
shapiro.test(df$faba)
hist(df$faba)
shapiro.test(df$pd)
hist(df$pd)
shapiro.test(df$gini)
hist(df$gini)
shapiro.test(df$v_s)
hist(df$v_s)


# Model initialization
model <- '
  # regressions
  produt ~ sr + ppt + faba + pd + gini + v_s
  '

# Model initialization
model <- '
 # regressions
 produt ~ sr + ppt + faba + pd + gini + v_s
 sr ~ ppt + v_s 
 gini ~ sr + ppt + v_s 
 pd ~ ppt + v_s + faba 
 '

# Correlation matrix, import ggcorrplot
install.packages("ggcorrplot")
library(ggcorrplot)
corr_matrix <- cor(df, use = "pairwise.complete.obs")

# Create the plot
ggcorr_plot <- ggcorrplot(corr_matrix, lab = TRUE)
ggcorr_plot + scale_fill_gradient2(low = "white", high = "black", mid = "grey", 
                                   midpoint = 0, limit = c(-1, 1), 
                                   space = "Lab", na.value = "grey50")


######## CHAT GPT #######3

library(lavaan)

# Carregar os dados
alldata <- readxl::read_excel("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/dadosmisto.xlsx")

# Selecionar as variáveis relevantes
variable_list <- c("log_produt", "SR", "ppt", "faba", "SESPD", "gini", "valor_s")
df <- alldata[variable_list]
colnames(df) <- c("produt", "sr", "ppt", "faba", "pd", "gini", "v_s")

# Especificação do modelo SEM
model <- ' 
  # Relações diretas
  produt ~ ppt + v_s
  
  # Relações indiretas
  sr ~ ppt + v_s
  gini ~ ppt + v_s + sr
  pd ~ ppt + v_s + faba
  produt ~ gini + pd + sr
  
  # Possíveis correlações baseadas nos dados e na teoria
  pd ~~ gini
  faba ~~ pd
'

# Ajustar o modelo
fit <- sem(model, data = df, estimator = "ML")

# Resumo dos resultados
summary(fit, standardized = TRUE, fit.measures = TRUE)

varTable(fit)

df_scaled <- scale(df)
fit <- sem(model, data = as.data.frame(df_scaled), estimator = "ML")
lavInspect(fit, "cov.lv")


# Modelo SEM mais simples: apenas efeitos diretos de ppt e v_s sobre produt
model_simple <- '
  # Relações diretas
  produt ~ ppt + v_s
'

# Ajuste do modelo
fit_simple <- sem(model_simple, data = df, estimator = "ML")

# Resumo dos resultados
summary(fit_simple, fit.measures = TRUE, standardized = TRUE)

## Modelo expandido
model_expanded <- '
  # Relações diretas
  produt ~ ppt + v_s + sr + gini + pd
  
  # Relações indiretas
  sr ~ ppt
  gini ~ ppt
  pd ~ ppt
'

# Ajuste do modelo expandido
fit_expanded <- sem(model_expanded, data = df, estimator = "ML")

# Resumo dos resultados
summary(fit_expanded, fit.measures = TRUE, standardized = TRUE)


## DEEP SEEK

# Verificar as variâncias
varTable(fit_expanded)

# Possivelmente padronizar as variáveis
df_std <- as.data.frame(scale(df))

model_revised <- '
  # Relações diretas (mantendo apenas significativas)
  produt ~ ppt + gini
  
  # Relações indiretas (mantendo apenas significativas/marginalmente significativas)
  sr ~ ppt
  gini ~ ppt
'

fit_revised <- sem(model_revised, data = df_std, estimator = "ML")
summary(fit_revised, fit.measures = TRUE, standardized = TRUE)

##

model_final <- '
  # Relações diretas
  produt ~ ppt + gini + sr  # incluir sr baseado na covariância residual
  
  # Mediações
  gini ~ ppt
  sr ~ ppt
  
  # Covariância entre termos de erro
  gini ~~ sr
'

fit_final <- sem(model_final, data = df_std, estimator = "ML")
summary(fit_final, fit.measures = TRUE, standardized = TRUE)

##

model_revised <- '
  # Relações diretas (remover sr não significativo)
  produt ~ ppt + gini
  
  # Mediações
  gini ~ ppt
  sr ~ ppt
  
  # Restrição: igualar efeitos de ppt
  ppt == a*ppt
'

fit_revised <- sem(model_revised, data = df_std)

##

model_better <- '
  # Relações diretas
  produt ~ ppt + gini
  
  # Mediações (mantendo apenas as significativas)
  gini ~ ppt
  
  # Deixamos sr como variável exógena correlacionada
  sr ~~ ppt + gini
'

fit_better <- sem(model_better, data = df_std, estimator = "ML")
summary(fit_better, standardized = TRUE)

##

model_mediation <- '
  # Variável dependente
  produt ~ c*ppt + b*gini
  
  # Mediador
  gini ~ a*ppt
  
  # Efeitos indiretos
  indirect := a*b
  total := c + (a*b)
'

fit_mediation <- sem(model_mediation, data = df_std)
summary(fit_mediation, standardized = TRUE)


##

# Efeitos padronizados totais
standardizedSolution(fit_better)

# Gráfico do modelo
# Versão mais robusta para gerar o diagrama
semPaths(fit_better, what = "std", 
         layout = "tree", 
         style = "lisrel",
         edge.label.cex = 0.8,
         sizeMan = 10,
         sizeLat = 10,
         nCharNodes = 4,
         rotation = 2,
         curve = 1.5,
         mar = c(3,3,3,3))

# Mostrar apenas caminhos significativos
semPaths(fit_better, "std", 
         thresholds = 0.05,  # só mostra p<0.05
         edge.label.cex = 0.8)

################# MODELO 2 ##################

alldata <- read_excel("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/dadosmisto.xlsx")

# Selecting relevant variables by their name in the dataset itself
variable_list <- c("log_produt", "SR", "ppt", "faba", "PSE", "gini","PC1nutri")
# Subsetting the dataframe
df <- alldata[variable_list]
# Renaming columns to make them easily referenced for SEM analysis
colnames(df) <- c("produt", "sr", "ppt", "faba", "PSE", "gini","nutri")

df_std <- as.data.frame(scale(df))

model_updated <- '
  # Relações diretas (variável dependente: produt)
  produt ~ ppt + gini + nutri
  
  # Mediações principais
  gini ~ ppt + nutri
  sr ~ ppt + nutri
  PSE ~ ppt + nutri + faba  # substituindo pd por PSE
  
  # Relações adicionais baseadas na teoria
  sr ~ gini
  PSE ~ sr
  
  # Covariâncias residuais
  gini ~~ sr
  PSE ~~ gini
'

# Ajuste o modelo com dados padronizados
fit_updated <- sem(model_updated, data = df_std, estimator = "ML")

# verifique se o modelo é teoricamente identificado: 
lavInspect(fit_updated, "free")$beta  # Matriz de coeficientes
lavInspect(fit_updated, "free")$psi   # Matriz de covariância

# Ajuste corrigido

model_identified <- '
  # Relações diretas (apenas 1 por variável dependente)
  produt ~ ppt + gini
  
  # Mediações (1 caminho por mediador)
  gini ~ nutri
  sr ~ ppt
  PSE ~ faba
  
  # Covariâncias restritas (apenas 1)
  gini ~~ PSE
  
  # Fixar variâncias residuais
  produt ~~ a*produt
  gini ~~ b*gini
  sr ~~ c*sr
  PSE ~~ d*PSE
  a > 0.001; b > 0.001; c > 0.001; d > 0.001
'

fit_identified <- sem(model_identified, 
                      data = df_std,
                      estimator = "ML",
                      se = "robust")
# 1. Verifique identificação
lavInspect(fit_identified, "info")$rank == length(lavInspect(fit_identified, "free"))

# 2. Resumo dos resultados
summary(fit_identified, standardized = TRUE, fit.measures = TRUE)

# 3. Verifique resíduos
residuals(fit_identified, type = "cor")$cov

## Soluções recomendadas

model_final_revised <- '
  # Efeitos diretos principais
  produt ~ ppt + gini
  
  # Apenas mediações significativas
  PSE ~ faba
  
  # Covariâncias significativas
  gini ~~ PSE
  sr ~~ PSE
  
  # Restrições de variância
  produt ~~ a*produt; gini ~~ b*gini
  sr ~~ c*sr; PSE ~~ d*PSE
  a > 0.001; b > 0.001; c > 0.001; d > 0.001
'

fit_final <- sem(model_final_revised, 
                 data = df_std,
                 estimator = "MLR")  # Usando erros robustos

# Verifique outliers
mahalanobis_dist <- lavInspect(fit_final, "residuals")$cov
which(mahalanobis_dist > qchisq(0.99, ncol(df_std)))

# Análise de poder
install.packages("semPower")
library(semPower)
semPower.aPriori(effect = 0.15, effect.measure = 'RMSEA', 
                 alpha = 0.05, power = 0.80, df = 5)

parameterEstimates(fit_final, 
                   boot.ci.type = "bca.simple") %>%
  filter(op == "~" | op == ":=")

library(semPlot)

# Diagrama básico garantido
semPaths(fit_final, what = "std", 
         layout = "tree",  # layout alternativo
         style = "lisrel",
         edge.label.cex = 0.7,
         sizeMan = 8,
         sizeLat = 8,
         nCharNodes = 4,
         rotation = 2)

semPaths(fit_final, "std",
         layout = "spring",
         groups = list(c("ppt","nutri","faba"), 
                       c("produt","gini","PSE")),
         pastel = TRUE)

##### PLS SEM #####

install.packages("semPLS")  # Instale se necessário
library(semPLS)

model_final_pls <- plsm(
  data = df_std,
  strucmod = list(
    # Relação principal de interesse
    c("produt", "PSE"),
    
    # Relações adicionais originais que você quer manter
    c("produt", "ppt"),
    c("produt", "gini"),
    
    # Covariâncias (se ainda relevantes)
    c("gini", "PSE")
  ),
  measuremod = list(
    produt = "produt",
    PSE = "PSE",
    ppt = "ppt",
    gini = "gini"
  )
)

# Estimar o modelo
fit_final <- plsr(model_final_pls, data = df_std)

# Resultados principais
results <- summary(fit_final)

# Coeficientes padronizados
path_coef <- results$path_coefficients
print(path_coef)

# Significância via bootstrap (500 amostras)
set.seed(123)
boot_results <- bootsempls(fit_final, nboot = 500)

### PLS SEM 2 ###

alldata <- read_excel("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/TRABALHO EXCEL/dadossem.xlsx")

# Selecting relevant variables by their name in the dataset itself
variable_list <- c("log_produt", "SR", "PC1_clima", "faba", "PSE", "gini","silte_g.kg")
# Subsetting the dataframe
df <- alldata[variable_list]
# Renaming columns to make them easily referenced for SEM analysis
colnames(df) <- c("produt", "sr", "aridez", "faba", "PSE", "gini","silte")

df_std <- as.data.frame(scale(df))

# --------------------------------------
# MODELO PLS-SEM - Efeitos diretos e indiretos
# Laíla Iglesias - Dados padronizados
# --------------------------------------

install.packages("plspm")
library(plspm)
library(readxl)

# 1. Carregar e preparar os dados
alldata <- read_excel("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/TRABALHO EXCEL/dadossem.xlsx")

# Selecionar e renomear variáveis
df <- alldata[c("log_produt", "SR", "PC1_clima", "faba", "PSE", "gini", "silte_g.kg")]
colnames(df) <- c("produt", "sr", "aridez", "faba", "pse", "gini", "silte")

# Padronizar os dados (média = 0, DP = 1)
df_std <- as.data.frame(scale(df))

# 2. Atualize a lista de blocos para usar os nomes EXATOS das colunas:
blocos <- list(
  aridez = "aridez",
  silte = "silte",
  faba = "faba",
  sr = "sr",
  gini = "gini",
  pse = "PSE",  # Note o uso de maiúsculas para corresponder ao dataframe
  produt = "produt"
)

# 3. Execute novamente o modelo:
resultado <- plspm(df_std, modelo, blocos, modes = modos,
                   boot.val = TRUE, br = 500)

# 6. Visualizar resultados
summary(resultado)

# 7. Plotar os coeficientes path
plot(resultado, what = "loadings")

# 8. Exportar resultados para CSV (opcional)
write.csv(resultado$path_coefs, "coeficientes_caminhos.csv")
write.csv(resultado$boot$paths, "significancia_caminhos.csv")

# Remove non-significant paths (e.g., aridez → produt)
modelo[7,1] <- 0  # Row 7 (produt), Column 1 (aridez)



modelo[4,1] <- 0  # Remove aridez → sr
modelo[6,1] <- 0  # Remove aridez → pse
modelo[7,6] <- 0  # Remove pse → produt

resultado <- plspm(df_std, modelo, blocos, modes = modos,
                   boot.val = TRUE, br = 500)
summary(resultado)


####### GRÁFICO #########

library(semPlot)
library(lavaan)

# 1. Create lavaan-style syntax from your significant paths
model_syntax <- '
# Structural model (only significant paths)
gini ~ a*aridez + b*silte + c*sr
pse ~ d*silte + e*faba
produt ~ f*silte + g*sr + h*gini

# Custom labels for coefficients
'

# 2. Fit a dummy lavaan model (we only need it for visualization)
fit <- sem(model_syntax, data=df_std)

# 3. Create the path diagram with proper formatting
semPaths(fit, what="est", 
         edge.label.cex=0.9,
         edge.color=ifelse(abs(coef(fit)) > 0.3, "red", "grey50"),
         layout="tree",
         sizeMan=8, 
         sizeLat=10,
         style="lisrel",
         nCharNodes=6,
         rotation=2,
         intercepts=FALSE,
         residuals=FALSE)

# Add title and legend
title("PLS-SEM Path Model\n(Red = |β| > 0.3, Grey = |β| ≤ 0.3)", line=1)
legend("bottomright", 
       legend=c("Strong Effect", "Weak Effect"), 
       col=c("red", "grey50"), 
       lty=1, lwd=2, bty="n")
