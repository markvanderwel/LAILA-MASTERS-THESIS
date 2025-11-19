

### DAVIES ET AL 2016 ###


library(ape)
library(picante)
library(lme4)
library(influence.ME)  # para DFBETA em LMM
library(ggtree)
library(dplyr)

install.packages("ggtree")

# 1) filogenia (você já criou tree_ok via V.PhyloMaker2)
library("V.PhyloMaker2")
library(writexl)

example <- read.csv("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/filogenia_total.csv")

tree <- phylo.maker(example, tree = GBOTB.extended.TPL,output.sp.list = TRUE,nodes = nodes.info.1.TPL, scenarios="S3")
tree_ok <- tree$scenario.3
summary(tree_ok)

# 2) dados por parcela (já carregados no seu script)
# dadosmisto1 <- ...

# 3) matriz comunidade (parcelas x espécies)
comunidade <- read.csv("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/comunidade_pd.csv",
                       row.names = 1, header = TRUE, sep = ";")

# manter apenas espécies que estão na árvore
sp_ok <- intersect(colnames(comunidade), tree_ok$tip.label)
comunidade <- comunidade[, sp_ok, drop = FALSE]

# (opcional) remover espécies da árvore que não estão na comunidade
tree_ok <- drop.tip(tree_ok, setdiff(tree_ok$tip.label, sp_ok))

### 2
library(picante)
library(tidyverse)

# Calcular SESPD (standardized effect size of PD)
ses_pd_vals <- ses.pd(samp = comunidade, tree = tree_ok,
                      null.model = "taxa.labels", runs = 999, include.root = FALSE)

# resultado tem as colunas: PD.obs, PD.rand.mean, PD.rand.sd, PD.obs.z, PD.obs.p, PD.obs.rank
# usamos o Z-score (PD.obs.z)

dados_pd <- dadosmisto |>
  mutate(SESPD1 = ses_pd_vals$pd.obs.z)

m_pd <- lmer(log_produt ~ SESPD1 + (1|site), data = dados_pd, REML = FALSE)
summary(m_pd)

## DFBETA

m_pd_lm <- lm(log_produt ~ SESPD1 + site, data = dados_pd)
dfb <- stats::dfbeta(m_pd_lm)[, "SESPD1"]   # nome exato da coluna

# Pacotes

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree")

library(ape)
library(phytools)   # getDescendants()
library(ggtree)
library(ggplot2)

# manter apenas espécies presentes na árvore
sp_ok <- intersect(colnames(comunidade), tree_ok$tip.label)
comunidade <- comunidade[, sp_ok, drop = FALSE]
tree_ok    <- drop.tip(tree_ok, setdiff(tree_ok$tip.label, sp_ok))

# (IMPORTANTE) mesma ordem de parcelas de dfb/comunidade/dados_pd
stopifnot(all(rownames(comunidade) == rownames(dados_pd)))
stopifnot(all(names(dfb) %in% rownames(comunidade)))


rownames(comunidade)[1:10]   # nomes das parcelas na matriz de comunidade
rownames(dados_pd)[1:10]     # nomes das parcelas nos dados de produtividade



dfb <- dfb[rownames(comunidade)]  # reordena o vetor de DFBETA

get_edge_incidence <- function(phy){
  E <- phy$edge
  nE <- nrow(E)
  ntip <- length(phy$tip.label)
  A <- matrix(0, nrow = ntip, ncol = nE,
              dimnames = list(phy$tip.label, paste0("e", seq_len(nE))))
  for(e in seq_len(nE)){
    nd <- E[e, 2]
    desc <- phytools::getDescendants(phy, nd)
    tips_idx <- intersect(desc, seq_len(ntip))
    if(length(tips_idx)){
      A[tips_idx, e] <- 1
    }
  }
  A
}

A_sp_edge <- get_edge_incidence(tree_ok)          # espécies × ramos (0/1)
PR <- (as.matrix(comunidade) %*% A_sp_edge) > 0   # parcelas × ramos (0/1)
PR <- PR * 1

# multiplica linha-a-linha cada parcela pelo seu DFBETA
PRw <- sweep(PR, 1, dfb, "*")                     # parcelas × ramos, ponderado
edge_influence <- colMeans(PRw, na.rm = TRUE)     # vetor por ramo


rownames(dados_pd) <- dados_pd$site
if (is.null(names(dfb))) names(dfb) <- rownames(dados_pd)

