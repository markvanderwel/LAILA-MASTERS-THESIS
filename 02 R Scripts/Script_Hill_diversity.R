# -----------------------------------------------
# ANÁLISE DE DIVERSIDADE FILOGENÉTICA BASEADA EM HILL NUMBERS
# -----------------------------------------------

# 1. CARREGAR PACOTES E DADOS
# Instalar pacotes se necessário (remova o # para instalar)
# install.packages(c("picante", "PhyloMeasures", "hillR", "ape", "phytools", "dplyr"))
install.packages("hillR")
library(picante)
library(hillR)
library(ape)
library(phytools)
library(dplyr)

# Carregar seus dados
# (Substitua pelos seus próprios dados)
comunidade <- read.csv("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/comunidade_pd.csv", row.names = 1,header = T, sep = ";")

# 2. CÁLCULO DAS MÉTRICAS

# a) Usando hillR (recomendado para simplicidade)
hill_metrics <- data.frame(
  Local = rownames(comunidade),
  PD_q0 = hill_phylo(comunidade, tree_ok, q = 0),
  PD_q1 = hill_phylo(comunidade, tree_ok, q = 1),
  PD_q2 = hill_phylo(comunidade, tree_ok, q = 2)
)

library(writexl)
setwd("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística")
write_xlsx(hill_metrics, "hill_metrics.xlsx")

# b) Função customizada para qPD(T) (Chao et al.)
calculate_qPD <- function(comunidade, tree_ok, q = 0, T = NULL) {
  if (is.null(T)) T <- max(node.depth.edgelength(tree_ok))
  
  branches <- tree_ok$edge
  branch_lengths <- tree_ok$edge.length
  node_ages <- node.depth.edgelength(tree_ok)
  root_age <- max(node_ages)
  
  valid_branches <- which((root_age - node_ages[branches[, 1]]) >= (root_age - T))
  
  node_abundances <- matrix(0, nrow = nrow(comm), ncol = max(branches))
  for (i in 1:nrow(comm)) {
    tip_abundances <- comm[i, ]
    for (tip in 1:length(tip_abundances)) {
      path <- nodepath(tree_ok, from = max(branches) + 1, to = tip)
      node_abundances[i, path[-length(path)]] <- node_abundances[i, path[-length(path)]] + tip_abundances[tip]
    }
  }
  
  qPD <- numeric(nrow(comm))
  for (i in 1:nrow(comm)) {
    a_i <- node_abundances[i, branches[valid_branches, 2]]
    L_i <- branch_lengths[valid_branches]
    
    if (q == 0) {
      qPD[i] <- sum(L_i)
    } else if (q == 1) {
      qPD[i] <- exp(-sum((L_i/T) * (a_i/sum(a_i)) * log(a_i/sum(a_i)))))
    } else {
      qPD[i] <- (sum(L_i * (a_i/sum(a_i))^q))^(1/(1-q))
    }
  }
  return(qPD)
}

# Calcular métricas customizadas
custom_metrics <- data.frame(
  Local = rownames(comunidade),
  qPD_0 = calculate_qPD(comunidade, tree_ok, q = 0),
  qPD_1 = calculate_qPD(comunidade, tree_ok, q = 1),
  qPD_2 = calculate_qPD(comunidade, tree_ok, q = 2)
)

# 3. Regressão

library(writexl)
library(readxl)
library(lmerTest) #MODELO MISTO
library(MuMIn) #AICc
setwd("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística")
dadosmisto <- read_excel("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/dadosmisto.xlsx")

# OFICIAL ### PRODUTIVIDADE ### site as random variable

# Modelo nulo sem os preditores fixos
modelo_nulo <- lmer(log_produt ~ (1 | site), data = dadosmisto, REML = FALSE)
summary(modelo_nulo)
AICc(modelo_nulo)


misto1 <- lmer(log_produt ~ PD_q0 + (1 | site), data = dadosmisto, REML = FALSE)
summary(misto1) # p-value = 0.0142 *
r.squaredGLMM(misto1) # 0.4845147
AICc(misto1) # 87.51476

misto2 <- lmer(log_produt ~ PD_q1 + (1 | site), data = dadosmisto, REML = FALSE)
summary(misto2) # p-value = 0.34
r.squaredGLMM(misto2) # 0.3855655
AICc(misto2) # 92.49965

misto3 <- lmer(log_produt ~ PD_q2 + (1 | site), data = dadosmisto, REML = FALSE)
summary(misto3) # p-value = 0.23
r.squaredGLMM(misto3) # 0.3933207
AICc(misto3) # 91.96249

misto4 <- lmer(log_produt ~ PD_q0 + SESPD + (1 | site), data = dadosmisto, REML = FALSE)
summary(misto4) # p-value = 0.23
r.squaredGLMM(misto4) # 0.3933207
AICc(misto4) # 91.96249


misto1b <- lmer(log_produt ~ poly(PD_q0, 2) + (1|site), dadosmisto)
summary(misto1b)

modelo_duplo <- lmer(log_produt ~ SR + PD_q0 + (1 | site), data = dadosmisto, REML = FALSE)
summary(modelo_duplo)

library(picante)

# 1. Verificar correspondência de nomes
nomes_arvore <- distinct$Species
nomes_comunidade <- colnames(comunidade)
setdiff(nomes_comunidade, nomes_arvore)
setdiff(nomes_arvore, nomes_comunidade)
all(nomes_comunidade %in% nomes_arvore)  # Deve ser TRUE

# Verificar se há comunidades vazias
comunidades_vazias <- rowSums(comunidade) == 0
sum(comunidades_vazias)  # Número de comunidades com abundância zero

# Verificar se há NA/NaN na matriz
any(is.na(comunidade))

# Checar correspondência exata (incluindo capitalização e caracteres especiais)
nomes_iguais <- all(colnames(comunidade) == distinct$Species[match(colnames(comunidade), distinct$Species)])

if (!nomes_iguais) {
  message("Atenção: Há diferenças nos nomes (capitalização, espaços, caracteres especiais)")
  
  # Mostrar as primeiras discrepâncias
  head(data.frame(
    Comunidade = colnames(comunidade),
    Arvore = distinct$Species[match(colnames(comunidade), distinct$Species)]
  ))
}

library(picante)

# Garantir que os nomes estão na mesma ordem
comunidade <- comunidade[, distinct$Species[distinct$Species %in% colnames(comunidade)]]

# Calcular distinct_medio com verificação explícita
dadosmisto$distinct_medio <- apply(comunidade, 1, function(x) {
  spp <- colnames(comunidade)[x > 0]
  pesos <- distinct$w[match(spp, distinct$Species)]
  
  # Verificação final
  if (any(is.na(pesos))) {
    stop("Espécies sem correspondência: ", paste(spp[is.na(pesos)], collapse = ", "))
  }
  
  weighted.mean(pesos, w = x[x > 0])
})
