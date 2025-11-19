## ============================
## 1) Pacotes
## ============================
# install.packages(c("ape","picante","phytools","ggplot2","dplyr","tibble","tidyr","scales"))
suppressPackageStartupMessages({
  library(ape)
  library(picante)
  library(phytools)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(scales)
})

## ggtree é do Bioconductor (opcional):
has_ggtree <- requireNamespace("ggtree", quietly = TRUE) && requireNamespace("treeio", quietly = TRUE)
if (has_ggtree) {
  library(ggtree)
}

## ============================
## 2) Leitura de dados
## ============================
# Ajuste para os seus caminhos se necessário:
# Filogenia foi gerada antes (tree_ok). Se já existir no ambiente, pule a leitura/criação.
# Exemplo: você já fez:
# example <- read.csv("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/filogenia_ml_regua.csv")
# tree <- phylo.maker(example, tree = GBOTB.extended.TPL, output.sp.list = TRUE, nodes = nodes.info.1.TPL, scenarios="S3")
# tree_ok <- tree$scenario.3

# Se já tem tree_ok no ambiente:
stopifnot(exists("tree_ok"))

# Traços:
dados <- read.csv("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/funcional_ml_regua.csv",
                  row.names = 1, header = TRUE, sep = ";", check.names = FALSE)

## ============================
## 3) Preparos (garantir binários, alinhamento e nomes)
## ============================
# Resolver eventuais nós politômicos (se existir):
tree_teste <- multi2di(tree_ok)  # transforma politomias em dicotomias mínimas

# Selecione os traços que quer analisar:
traits_raw <- dados[, c("WD","SLA","LDMC")]
# garantir numéricos:
traits_raw <- mutate_all(traits_raw, ~ suppressWarnings(as.numeric(.)))

# Remover linhas com NA em TODOS os traços (ou troque por rowSums(is.na(.)) < ncol(.)):
traits_raw <- traits_raw[complete.cases(traits_raw), , drop = FALSE]

# Interseção de espécies:
sp_common <- intersect(rownames(traits_raw), tree_teste$tip.label)
if (length(sp_common) < 3) stop("Poucas espécies em comum entre traços e árvore.")

# Reduzir árvore e matriz:
tree_use   <- drop.tip(tree_teste, setdiff(tree_teste$tip.label, sp_common))
traits_use <- traits_raw[sp_common, , drop = FALSE]

# Ordenar a matriz pela ordem dos ápices:
traits_use <- traits_use[tree_use$tip.label, , drop = FALSE]

## ============================
## 4) Função robusta para Blomberg's K
## ============================
get_sig <- function(traits, tree) {
  tr_names <- colnames(traits)
  res <- lapply(tr_names, function(tr) {
    x <- setNames(traits[, tr], rownames(traits))
    x <- x[tree$tip.label]  # garante ordem
    out <- try(phytools::phylosig(tree, x, method = "K", test = TRUE), silent = TRUE)
    if (inherits(out, "try-error") || is.null(out$K)) {
      data.frame(trait = tr, k = NA_real_, p = NA_real_)
    } else {
      data.frame(trait = tr, k = out$K, p = out$P)
    }
  })
  dplyr::bind_rows(res)
}

tbl_sig <- get_sig(traits_use, tree_use)
print(tbl_sig)

## ============================
## 5) Figura A — Barras de K
## ============================
gg_sig <- ggplot(tbl_sig, aes(x = reorder(trait, k), y = k)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_text(aes(label = ifelse(!is.na(p) & p < 0.05, "*", "")),
            vjust = -0.5, size = 6) +
  coord_cartesian(ylim = c(0, max(1.05, max(tbl_sig$k, na.rm = TRUE) * 1.15))) +
  labs(x = NULL, y = "Blomberg's K",
       title = "Phylogenetic Signal by Trait (K)",
       subtitle = "* denotes p < 0.05 (randomization test)") +
  theme_minimal(base_size = 12)
print(gg_sig)

setwd("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/Resultados funcional/Gráficos")
ggsave("fig_K_por_traco.png", gg_sig, width = 5, height = 4, dpi = 300)

## ============================
## 6) Figura B — contMap (um por traço)
## ============================
# Dica: salve cada mapa se quiser inserir no slide
for (tr in colnames(traits_use)) {
  x <- setNames(traits_use[, tr], rownames(traits_use))[tree_use$tip.label]
  cm <- contMap(tree_use, x, plot = FALSE)
  # Abrir dispositivo gráfico de alta resolução (opcional):
  png(paste0("contmap_", tr, ".png"), width = 1600, height = 1200, res = 200)
  plot(cm, fsize = 0.5, legend = 0.7*max(nodeHeights(tree_use)), outline = FALSE)
  title(tr)
  dev.off()
}

## ============================
## 7) Figura C — Heatmap de traços nos ápices
## ============================
## ===== Figura C — Heatmap de traços nos ápices =====
# Padroniza (z-score) e garante ordem dos ápices
traits_scaled <- scale(traits_use) |> as.data.frame()
traits_scaled$species <- rownames(traits_use)
traits_scaled <- traits_scaled[match(tree_use$tip.label, traits_scaled$species), ]
rownames(traits_scaled) <- traits_scaled$species

if (requireNamespace("ggtree", quietly = TRUE) && requireNamespace("treeio", quietly = TRUE)) {
  ggtree::ggtree(tree_use) + theme_tree2() -> p_tree
  p_heat <- ggtree::gheatmap(p_tree,
                             traits_scaled[, c("WD","SLA","LDMC")],
                             offset = 0.02, width = 0.35, colnames = TRUE) +
    scale_fill_gradient2(low = scales::muted("blue"), mid = "white", high = scales::muted("red"),
                         midpoint = 0, name = "z-score") +
    theme(legend.position = "right")
  print(p_heat)
  # ggsave("heatmap_traits_tip.png", p_heat, width = 9, height = 7, dpi = 300)
  
} else {
  # Alternativa SEM ggtree: dendrograma + heatmap
  hc <- as.hclust.phylo(tree_use)
  dend <- as.dendrogram(hc)
  
  long_df <- traits_scaled |>
    dplyr::select(species, WD, SLA, LDMC) |>
    tidyr::pivot_longer(cols = c(WD, SLA, LDMC), names_to = "trait", values_to = "z")
  
  long_df$species <- factor(long_df$species, levels = tree_use$tip.label)
  
  # Dendrograma (base R)
  plot(dend, horiz = TRUE, main = "Filogenia (dendrograma)")
  
  # Heatmap (ggplot2)
  gg_heat <- ggplot(long_df, aes(x = trait, y = species, fill = z)) +
    geom_tile() +
    scale_fill_gradient2(low = scales::muted("blue"), mid = "white", high = scales::muted("red"),
                         midpoint = 0, name = "z-score") +
    labs(x = NULL, y = NULL, title = "Heatmap de traços nos ápices (z-score)") +
    theme_minimal(base_size = 12) +
    theme(axis.text.y = element_text(size = 6))
  print(gg_heat)
  # ggsave("heatmap_traits_tip_no_ggtree.png", gg_heat, width = 6, height = 10, dpi = 300)
}


## ============================
## 8) Saída resumida no console
## ============================
cat("\n=== Resumo K e p ===\n")
print(tbl_sig %>% arrange(desc(k)))
cat("\nInterpretação rápida: K > 1 sugere sinal filogenético mais forte que o esperado sob BM; * indica p < 0.05.\n")
