
# ---- SEM ----

library(lavaan)

# Load the data

dadosmisto <- read.csv("01 Datasets/01_raw_data/dadosmisto.csv",
                       header = TRUE)

dadosmisto1 <- dadosmisto[-c(1:7), ]

# Scale

dados_scaled <- as.data.frame(scale(dadosmisto1[, sapply(dadosmisto1, is.numeric)]))

# ---- Model A ----

model_A <- '
  ############################
  # Climate → Soil
  ############################
  c.n_soloid ~ season_ppt

  ############################
  # Phylogenetic/Functional structure
  ############################
  sr    ~ faba
  pcps1 ~ pse + faba

  ############################
  # Soil + phylo → Biomass
  ############################
  log_produt ~ sr + pcps1 + c.n_soloid + season_ppt

  ############################
  # Correlated exogenous
  ############################
  faba ~~ pse
  season_ppt ~~ faba + pse
'


fit_finalA <- sem(
  model = model_A,
  data  = dados_scaled,
  estimator = "ML"
)

summary(fit_finalA, standardized = TRUE, fit.measures = TRUE)

# ---- Model B ----

model_B <- '
  ############################
  # Climate → Soil
  ############################
  c.n_soloid ~ season_ppt

  ############################
  # Soil influences structure
  ############################
  sr    ~ faba + c.n_soloid
  pcps1 ~ pse + faba + c.n_soloid

  ############################
  # Biomass drivers
  ############################
  log_produt ~ sr + pcps1 + c.n_soloid + season_ppt

  ############################
  # Correlated exogenous
  ############################
  faba ~~ pse
  season_ppt ~~ faba + pse
'

fit_finalB <- sem(
  model = model_B,
  data  = dados_scaled,
  estimator = "ML"
)

summary(fit_finalB, standardized = TRUE, fit.measures = TRUE)

# ---- Model C ----

model_C <- '
  ############################
  # Climate affects phylogenetic structure
  ############################
  sr    ~ faba + season_ppt
  pcps1 ~ pse + faba + season_ppt

  ############################
  # Climate → Soil
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
'
fit_finalC <- sem(
  model = model_C,
  data  = dados_scaled,
  estimator = "ML"
)

summary(fit_finalC, standardized = TRUE, fit.measures = TRUE)
