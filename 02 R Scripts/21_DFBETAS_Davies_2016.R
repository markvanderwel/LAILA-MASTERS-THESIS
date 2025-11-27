library(readr)
library(tidyverse)
library("V.PhyloMaker2")
library(picante)

com_matrix <- read.csv("01 Datasets/01_raw_data/comunidade_abund.csv",
                       row.names = 1,
                       sep = ";")

com_matrix <- read.csv(
  "01 Datasets/01_raw_data/comunidade_abund.csv",
  row.names = 1,
  sep = ";",
  check.names = FALSE
)
# forçar tudo para numeric (mantendo nomes)
com_matrix <- as.data.frame(
  lapply(com_matrix, function(x) as.numeric(as.character(x))),
  row.names = rownames(com_matrix)
)

dadosmisto <- read.csv("01 Datasets/01_raw_data/dadosmisto.csv",
                       header = TRUE)


## Phylo tree

example <- read_csv("~/01 Masters_LA/00 MASTERS-DATA/01 Datasets/01_raw_data/filogenia_total.csv")

tree <- phylo.maker(example, tree = GBOTB.extended.TPL,output.sp.list = TRUE,nodes = nodes.info.1.TPL, scenarios="S3")
tree_ok <- tree$scenario.3

###


m4 <- lmer(log_produt ~ pd + (1 | site), data = dadosmisto, REML = FALSE)
summary(m4)

install.packages("influence.ME")
library(influence.ME)


# Calcular influência por observação
infl <- influence(m4, obs = TRUE)


dfbeta_vals <- dfbetas(infl)
head(dfbeta_vals)

dfbeta_df <- as.data.frame(dfbeta_vals)
dfbeta_df$parcela <- dadosmisto$parcela
head(dfbeta_df)


dfbeta_pd <- dfbeta_df$pd
names(dfbeta_pd) <- dfbeta_df$parcela

# 1. Ver quais parcelas existem nos dois objetos
common_plots <- intersect(rownames(com_matrix), names(dfbeta_pd))

length(common_plots)  # só para conferir

# 2. Restringir comunidade e DFBETA a esse conjunto comum
com_matrix_use  <- com_matrix[common_plots, , drop = FALSE]
dfbeta_pd_use   <- dfbeta_pd[common_plots]

# 3. Conferir se a ordem bate
identical(rownames(com_matrix_use), names(dfbeta_pd_use))

# Peso por parcela = |DFBETA de pd|
w_plot <- abs(dfbeta_pd_use[rownames(com_matrix_use)])

# Matriz de presença/ausência
pres_mat <- com_matrix_use > 0

# Peso de cada espécie = soma dos pesos das parcelas onde ela ocorre
spp_weights <- colSums(pres_mat * w_plot)

summary(spp_weights)
range(spp_weights)

install.packages("phangorn")   # se ainda não tiver
library(phangorn)

# Pesos na ordem de tree_ok$tip.label
tip_w <- spp_weights[match(tree_ok$tip.label, names(spp_weights))]
tip_w[is.na(tip_w)] <- 0  # espécies da árvore que não aparecem em nenhuma parcela

branch_weights <- numeric(nrow(tree_ok$edge))

for (e in seq_len(nrow(tree_ok$edge))) {
  child_node <- tree_ok$edge[e, 2]
  
  # todos os tips descendentes desse nó
  desc_tips <- Descendants(tree_ok, child_node, type = "tips")[[1]]
  
  # soma os pesos das espécies descendentes
  branch_weights[e] <- sum(tip_w[desc_tips])
}

names(branch_weights) <- apply(tree_ok$edge, 1, paste, collapse = "-")

summary(branch_weights)
range(branch_weights)

# Normalizar para 0–1
bw_norm <- branch_weights / max(branch_weights)
bw_norm[is.na(bw_norm)] <- 0

pal <- colorRampPalette(c("lightgrey", "black"))
cols <- pal(100)
edge_cols <- cols[1 + floor(99 * bw_norm)]

plot(tree_ok, edge.color = edge_cols, show.tip.label = FALSE)

library(dplyr)
library(phangorn)

# ordenar galhos do maior pro menor peso
ord_edges <- order(branch_weights, decreasing = TRUE)

# ver um resuminho dos 10 mais importantes
branch_weights[ord_edges[1:10]]

top_n <- 10  # número de galhos que você quer inspecionar

lista_clados <- lapply(ord_edges[1:top_n], function(e) {
  child_node <- tree_ok$edge[e, 2]
  tips <- Descendants(tree_ok, child_node, type = "tips")[[1]]
  
  data.frame(
    edge_id = e,
    branch_weight = branch_weights[e],
    species = tree_ok$tip.label[tips],
    stringsAsFactors = FALSE
  )
})

clados_top <- bind_rows(lista_clados, .id = "clade_rank")

View(clados_top)

jpeg("~/01 Masters_LA/06 Figures/03 Plot_phylo_trees/pd_dfbeta_branchweights.jpeg",
     width = 3000, height = 3000, res = 300)

par(mar = c(1,1,1,1))
plot(tree_ok, edge.color = edge_cols, show.tip.label = TRUE, cex = 0.4)

dev.off()


fabaceae <- c(
  "Anadenanthera_colubrina",
  "Andira_anthelmia",
  "Apuleia_leiocarpa",
  "Bauhinia_forficata",
  "Cenostigma_pluviosum",
  "Centrolobium_microchaete",
  "Centrolobium_tomentosum",
  "Chloroleucon_tortum",
  "Copaifera_langsdorffii",
  "Dalbergia_nigra",
  "Enterolobium_contortisiliquum",
  "Erythrina_verna",
  "Hydrochorea_pedicellaris",
  "Hymenaea_courbaril",
  "Inga_edulis",
  "Inga_laurina",
  "Inga_vera",
  "Libidibia_ferrea",
  "Machaerium_nyctitans",
  "Marlimorimia_contorta",
  "Mimosa_bimucronata",
  "Mimosa_schomburgkii",
  "Paubrasilia_echinata",
  "Peltophorum_dubium",
  "Pityrocarpa_gonoacantha",
  "Pityrocarpa_paniculata",
  "Plathymenia_reticulata",
  "Pterocarpus_rohrii",
  "Pterogyne_nitens",
  "Schizolobium_parahyba",
  "Senegalia_polyphylla",
  "Senna_macranthera",
  "Senna_multijuga",
  "Stryphnodendron_adstringens"
)

library(phangorn)

# índice do galho mais pesado
ord_edges  <- order(branch_weights, decreasing = TRUE)
edge_sel   <- ord_edges[1]
child_node <- tree_ok$edge[edge_sel, 2]

# tips descendentes desse galho
tips_sel <- Descendants(tree_ok, child_node, type = "tips")[[1]]
species_clade <- tree_ok$tip.label[tips_sel]

species_clade

# quantas das espécies do clado são Fabaceae?
sum(species_clade %in% fabaceae)
length(species_clade)
prop_fab <- sum(species_clade %in% fabaceae) / length(species_clade)

prop_fab
intersect(species_clade, fabaceae)


# Conclusions

# Although PD was significant in the models, this effect disappeared when PD was standardized by species richness (SES-PD), indicating that PD was acting largely as a proxy for richness. Thus, species richness itself is an important driver of biomass accumulation in these restoration plantings.However, PCPS revealed that certain phylogenetic lineages—particularly Fabaceae—also contributed disproportionately to biomass production. Therefore, biomass in restored forests is influenced both by higher species richness and by the presence of specific functionally important clades. These results suggest that the most effective restoration strategy is not only increasing species richness, but also ensuring representation of key evolutionary lineages, such as nitrogen-fixing Fabaceae.”