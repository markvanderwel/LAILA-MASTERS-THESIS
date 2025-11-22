
############################################
# Phylogenetic signal (K) and trait maps
# Sites: ML + REGUA
# Author: Laila Iglesias
############################################

## 1) Load packages ----
library(ape)
library(phytools)
library(dplyr)
library(ggplot2)
library(V.PhyloMaker2)

## 2) Read species list and traits ----

# Species list for the phylogeny
example <- read.csv(
  "~/01 Masters_LA/00 MASTERS-DATA/01 Datasets/01_raw_data/filogenia_ml_regua.csv"
)

# Traits (rows = species, columns = traits)
dados <- read.csv(
  "~/01 Masters_LA/00 MASTERS-DATA/01 Datasets/01_raw_data/funcional_ml_regua.csv",
  row.names = 1,
  header = TRUE,
  sep = ";"
)

## 3) Build phylogenetic tree ----

tree_full <- phylo.maker(
  sp.list = example,
  tree = GBOTB.extended.TPL,
  nodes = nodes.info.1.TPL,
  scenarios = "S3"
)

tree_ok <- tree_full$scenario.3
tree_ok <- multi2di(tree_ok)  # resolve polytomies

## 4) Prepare trait matrix ----

# Select traits of interest
traits_raw <- dados[, c("WD", "SLA", "LDMC")]

# Make sure everything is numeric
traits_raw <- mutate_all(traits_raw, as.numeric)

# Keep only species that are both in the tree and in the trait table
sp_common <- intersect(rownames(traits_raw), tree_ok$tip.label)

tree_use   <- drop.tip(tree_ok, setdiff(tree_ok$tip.label, sp_common))
traits_use <- traits_raw[sp_common, , drop = FALSE]

# Order traits by tree tip order
traits_use <- traits_use[tree_use$tip.label, , drop = FALSE]

## 5) Blomberg's K for each trait ----

res_k <- data.frame(
  trait = colnames(traits_use),
  K     = NA_real_,
  p     = NA_real_
)

for (i in seq_along(res_k$trait)) {
  tr_name <- res_k$trait[i]
  
  x <- traits_use[, tr_name]
  names(x) <- rownames(traits_use)
  x <- x[tree_use$tip.label]
  
  k_out <- phylosig(tree_use, x, method = "K", test = TRUE)
  
  res_k$K[i] <- k_out$K
  res_k$p[i] <- k_out$P
}

print(res_k)

## 6) Plot – barplot of K ----

fig_k <- ggplot(res_k, aes(x = trait, y = K)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_text(
    aes(label = ifelse(p < 0.05, "*", "")),
    vjust = -0.5,
    size = 5
  ) +
  labs(
    x = NULL,
    y = "Blomberg's K",
    title = "Phylogenetic signal by trait",
    subtitle = "* indicates p < 0.05 (randomization test)"
  ) +
  theme_minimal(base_size = 12)

print(fig_k)

ggsave(
  "~/01 Masters_LA/06 Figures/04 Plots_Functional_Diversity/fig_K_por_traco.png",
  fig_k,
  width = 5,
  height = 4,
  dpi = 300
)

## 7) Trait maps on the phylogeny (contMap) ----

output_dir <- "~/01 Masters_LA/06 Figures/04 Plots_Functional_Diversity/"

for (tr in colnames(traits_use)) {
  x <- traits_use[, tr]
  names(x) <- rownames(traits_use)
  x <- x[tree_use$tip.label]
  
  cm <- contMap(tree_use, x, plot = FALSE)
  
  png(
    filename = paste0(output_dir, "contmap_", tr, ".png"),
    width = 1600,
    height = 1200,
    res = 200
  )
  
  plot(
    cm,
    fsize  = 0.5,
    legend = 0.7 * max(nodeHeights(tree_use)),
    outline = FALSE
  )
  title(tr)
  
  dev.off()
}

