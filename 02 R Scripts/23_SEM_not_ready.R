#############
#### SEM ####
#############


# Install and load lavaan if not already installed
install.packages("lavaan")
library(lavaan)

# Load the data

dadosmisto <- read.csv("01 Datasets/01_raw_data/dadosmisto.csv",
                       header = TRUE)

dadosmisto1 <- dadosmisto[-c(1:7), ]

# Define the SEM model
model <- '
  season_ppt ~ altitude
  c.n_soloid ~ season_ppt
  sr ~ season_ppt
  pcps1 ~ season_ppt
  log_produt ~ c.n_soloid + sr + pcps1 + altitude + season_ppt
'


# Fit the model (assuming your data frame is called df)
fit <- sem(model, data = dadosmisto1)

# Summary of results
summary(fit, standardized = TRUE, fit.measures = TRUE)

######################################

model_revised <- '
  # Keep essential direct predictors of log_produt (based on lmer m12 and m17 results)
  # These are the *only* direct effects on log_produt
  log_produt ~ sr + pcps1 + c.n_soloid + season_ppt

  # Define potential mediators: How do these predictors relate to upstream variables?
  # The lmer models tested direct effects, the SEM now models the *system*.

  # Keep the strong, significant upstream links from the original model:
  c.n_soloid ~ season_ppt  # Strong link in original SEM

  # The correlation matrix image suggests many other possible links
  # (e.g., season_ppt is correlated with altitude, Tmin, Tmax, etc.)
  # We should include only theoretically relevant mediators based on the lmer findings.
'

# Fit the model (assuming your data frame is called df)
fit <- sem(model_revised, data = dadosmisto1)

summary(fit, standardized = TRUE, fit.measures = TRUE)


###################################################

model_identified <- '
  # Regressions (Causal paths remain the same)
  log_produt ~ sr + pcps1 + c.n_soloid + season_ppt
  c.n_soloid ~ season_ppt

  # Covariances (Correlations only between truly independent variables)
  # SR and PCPS1 are exogenous predictors
  sr ~~ pcps1
  
  # SR and Season_PPT are exogenous predictors
  sr ~~ season_ppt
  
  # PCPS1 and Season_PPT are exogenous predictors
  pcps1 ~~ season_ppt
  
  # Note: You do *not* add c.n_soloid ~~ season_ppt here,
  # because you already specified a causal path between them above.
'

fit_identified <- sem(model_identified, data = dadosmisto1)
summary(fit_identified, standardized = TRUE, fit.measures = TRUE)

### ONLY DIVERSITY ###

library(lavaan)

# Path model: Fabaceae → SR → Biomass
#             PSE → PCPS1 → Biomass
model_phylo <- '
  # Regressions (setas do diagrama)
  sr        ~ faba         # Fabaceae → Species richness
  pcps1     ~ pse          # PSE → PCPS1
  log_produt ~ sr + pcps1  # SR and PCPS1 → Biomass accumulation

  # (Optional) Correlation between Fabaceae and PSE
  # Uncomment if you want to allow them to covary:
  # faba ~~ pse
'

fit_phylo <- sem(
  model = model_phylo,
  data  = dadosmisto,
  estimator = "ML"
)

summary(fit_phylo, standardized = TRUE, fit.measures = TRUE)

model_phylo_direct <- '
  # Regressions
  sr        ~ faba
  pcps1     ~ pse
  log_produt ~ sr + pcps1 + faba + pse

  # Optional covariance
  # faba ~~ pse
'

fit_phylo_direct <- sem(
  model = model_phylo_direct,
  data  = dadosmisto,
  estimator = "ML"
)

summary(fit_phylo_direct, standardized = TRUE, fit.measures = TRUE)

#####################################################


library(lavaan)

model_phylo_env <- '
  ############################
  # 1) Phylogenetic structure
  ############################
  sr    ~ faba        # Fabaceae → richness
  pcps1 ~ pse         # PSE → PCPS1

  ############################
  # 2) Soil mediated by climate
  ############################
  c.n_soloid ~ season_ppt     # rainfall seasonality → soil C:N

  ############################
  # 3) Biomass accumulation
  ############################
  log_produt ~ sr + pcps1 + c.n_soloid + season_ppt

  ############################
  # 4) Correlated exogenous variables
  ############################
  faba ~~ pse
  season_ppt ~~ faba + pse
'

fit_phylo_env <- sem(
  model = model_phylo_env,
  data  = dadosmisto,
  estimator = "ML"
)

summary(fit_phylo_env, standardized = TRUE, fit.measures = TRUE)

### STANDARDIZE

library(lavaan)

###############################
# 1. Scale das variáveis
###############################

vars_to_scale <- c(
  "sr","pcps1","pse","faba","c.n_soloid",
  "season_ppt","log_produt"
)

dados_scaled <- as.data.frame(scale(dadosmisto[, vars_to_scale]))


###############################
# 2. Modelo final reduzido
###############################

model_phylo_env_scaled <- '
  ############################
  # Phylogenetic structure
  ############################
  sr    ~ faba        # Fabaceae → richness
  pcps1 ~ pse         # PSE → PCPS1

  ############################
  # Soil mediated by climate
  ############################
  c.n_soloid ~ season_ppt     # rainfall seasonality → soil C:N

  ############################
  # Biomass accumulation
  ############################
  log_produt ~ sr + pcps1 + c.n_soloid + season_ppt

  ############################
  # Correlated exogenous variables
  ############################
  faba ~~ pse
  season_ppt ~~ faba + pse
'

###############################
# 3. Rodar o SEM
###############################

fit_phylo_env_scaled <- sem(
  model = model_phylo_env_scaled,
  data  = dados_scaled,
  estimator = "ML"
)

###############################
# 4. Resultados
###############################

summary(fit_phylo_env_scaled,
        standardized = TRUE,
        fit.measures = TRUE)

modindices(fit_phylo_env_scaled, sort = TRUE)[1:20, ]


################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

model_final <- '
  ############################
  # Phylogenetic structure
  ############################
  sr    ~ faba        # Fabaceae → richness
  pcps1 ~ pse + faba  # PSE + Fabaceae → PCPS1 (ESSENCIAL)

  ############################
  # Soil mediated by climate
  ############################
  c.n_soloid ~ season_ppt

  ############################
  # Biomass accumulation
  ############################
  log_produt ~ sr + pcps1 + c.n_soloid + season_ppt

  ############################
  # Correlated exogenous
  ############################
  faba ~~ pse
  season_ppt ~~ faba + pse
'

fit_final <- sem(
  model = model_final,
  data  = dados_scaled,
  estimator = "ML"
)

summary(fit_final, standardized = TRUE, fit.measures = TRUE)


# Plot
install.packages("semPlot")
library(semPlot)


# Instalar pacotes se necessário
# install.packages("lavaan")
# install.packages("semPlot")


library(lavaan)
library(semPlot)

# Supondo que você já tenha o modelo ajustado: fit_final

# Nomes bonitos
labels_bonitos <- c(
  "faba"        = "Fator ABA",
  "pse"         = "PSE",
  "season_ppt"  = "Precipitação Sazonal",
  "sr"          = "Riqueza de Espécies",
  "pcps1"       = "PCPS1",
  "c.n_soloid"  = "Carbono no Solo",
  "log_produt"  = "Produtividade (log)"
)

# Gerar diagrama SEM
semPaths(
  object = fit_final,
  what = "std",                # coeficientes padronizados
  layout = "tree",             # layout em árvore
  style = "lisrel",            # estilo clássico SEM
  edge.label.cex = 1.0,        # tamanho dos coeficientes
  sizeMan = 8,                 # tamanho das caixas
  sizeLat = 10,
  nCharNodes = 0,              # não cortar nomes
  nodeLabels = labels_bonitos, # aplicar nomes bonitos
  residuals = TRUE,            # mostrar resíduos
  intercepts = FALSE,          # não mostrar interceptos
  edge.color = "black",        # cor padrão (vamos usar edge.color = "black")
  color = list(lat = "white", man = "lightgray"), # cores das caixas
  edge.width = 2,              # espessura base
  mar = c(5, 5, 5, 5),
  curvePivot = TRUE
)

# Para cores automáticas e espessura proporcional:
# Use argument 'edge.color' = "std" e 'edge.width' = "std"
# Exemplo:
semPaths(
  fit_final,
  what = "std",
  layout = "tree",
  style = "lisrel",
  edge.label.cex = 1.0,
  sizeMan = 8,
  nodeLabels = labels_bonitos,
  residuals = TRUE,
  intercepts = FALSE,
  edge.color = "std",   # cores baseadas no sinal (verde/vermelho)
  edge.width = "std",   # espessura proporcional ao coeficiente
  color = list(lat = "white", man = "lightgray"),
  mar = c(5, 5, 5, 5)
)
