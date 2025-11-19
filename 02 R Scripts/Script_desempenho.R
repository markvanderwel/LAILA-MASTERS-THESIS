################ DESEMPENHO POR ESPÉCIE #########################

# 1. Load packages
library(readxl)
library(tidyverse)
library(writexl)
library(lmerTest)  

# 2. Carregar os dados

dados <- read_excel("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/dados_desempenho_semcsa.xlsx")

lmer(produt_z_g.ano ~ especie * PC1_nutri + (1|site.parcela), data = dados)

library(ggplot2)
library(dplyr)

# Escolher um subconjunto de espécies para não poluir o gráfico (ex: 8 mais frequentes)

top_species <- dados %>%
  count(especie, sort = TRUE) %>%
  slice_max(order_by = n, n = 8) %>%
  pull(especie)


dados_subset <- dados %>%
  filter(especie %in% top_species)

# Plot
ggplot(dados_subset,
       aes(x = PC1_nutri, y = produt_z_g.ano, color = especie)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_y_continuous(trans = "log10") +   # << log to “decompress”
  facet_wrap(~ especie, scales = "free") + # free x and y scales per panel
  labs(title = "Species responses to the climatic gradient (PC1)",
       x = "Climatic gradient (PC1_clima)",
       y = "Biomass (z-score, log10 scale)") +
  theme_minimal()

# Set working directory
setwd("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/graficos")

# Save with ggsave
ggsave("species_response.jpeg",
       width = 2000, height = 1600, units = "px", dpi = 300)

##### forest plot of species-specific slopes (productivity vs. PC1_clima #####

# ==== PACKAGES ====
library(readxl)
library(dplyr)
library(purrr)
library(broom)
library(ggplot2)
library(forcats)
library(stringr)

# ==== INPUTS ====
# Adjust paths and sheet name as needed:
excel_path <- "C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/dados_desempenho.xlsx"
sheet_name <- "principal"
out_dir    <- "C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/graficos"
min_n      <- 10   # increase to 20 for more stable CIs

# Column name for productivity (use your exact column name)
prod_col <- "produt_z_g.ano"

# ==== LOAD DATA ====
dat <- read_xlsx(excel_path, sheet = sheet_name) |>
  mutate(
    especie   = as.character(especie),
    PC1_clima = as.numeric(PC1_clima)
  )

stopifnot(prod_col %in% names(dat))

# ==== HELPER: fit slopes per species (simple lm: y ~ x) ====
fit_slopes <- function(df, response, predictor = "PC1_clima", min_n = 10) {
  df |>
    select(especie, !!sym(predictor), !!sym(response)) |>
    rename(x = !!sym(predictor), y = !!sym(response)) |>
    filter(is.finite(x), is.finite(y)) |>
    group_by(especie) |>
    filter(n() >= min_n, var(x) > 0) |>      # ensure enough data & variation in x
    group_modify(~{
      m <- lm(y ~ x, data = .x)
      tb <- tidy(m)
      tibble(
        n       = nrow(.x),
        slope   = coef(m)[["x"]],
        lwr     = confint(m)["x", 1],
        upr     = confint(m)["x", 2],
        p_value = tb$p.value[tb$term == "x"]
      )
    }) |>
    ungroup() |>
    arrange(slope) |>
    mutate(especie = fct_inorder(especie))
}

# ==== HELPER: forest plot ====
plot_forest <- function(slopes_tbl, title, xlab) {
  ggplot(slopes_tbl,
         aes(x = slope, y = fct_reorder(especie, slope))) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2) +
    geom_point(size = 2) +
    labs(title = title,
         x = xlab,
         y = "Species") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white"),
      plot.background  = element_rect(fill = "white"),
      axis.text        = element_text(color = "black"),
      axis.title       = element_text(color = "black"),
      title            = element_text(color = "black")
    )
}

# ==== PRODUCTIVITY ====
slopes_prod <- fit_slopes(dat, response = prod_col, min_n = min_n)

forest_prod <- plot_forest(
  slopes_prod,
  title = "Species-specific responses to aridity (PC1): Productivity",
  xlab  = "Slope of productivity (z-score) vs PC1_clima (β)"
)

# ==== SAVE OUTPUTS ====
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(out_dir, "forest_slopes_PC1_productivity.jpeg"),
       forest_prod,
       width = 8,
       height = max(6, 0.4 * nrow(slopes_prod)),
       units = "in",
       dpi = 300)

write.csv(slopes_prod,
          file.path(out_dir, "forest_slopes_PC1_productivity.csv"),
          row.names = FALSE)

# ==== OPTIONAL: quick variant with stricter sample size ====
# slopes_prod20 <- fit_slopes(dat, response = prod_col, min_n = 20)
# forest_prod20 <- plot_forest(
#   slopes_prod20,
#   title = "Responses to aridity (PC1): Productivity (n ≥ 20)",
#   xlab  = "Slope of productivity (z-score) vs PC1_clima (β)"
# )
# ggsave(file.path(out_dir, "forest_slopes_PC1_productivity_n20.png"),
#        forest_prod20, width = 8,
#        height = max(6, 0.4 * nrow(slopes_prod20)),
#        units = "in", dpi = 300)


###################### Nutrients ###############################

# 1. Load packages
library(readxl)
library(tidyverse)
library(writexl)
library(lmerTest)  
library(ggplot2)
library(dplyr)

# 2. Carregar os dados

dados <- read_excel("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/dados_desempenho.xlsx")

lmer(produt_z_g.ano ~ especie * PC1_nutri + (1|site.parcela), data = dados)

# Escolher um subconjunto de espécies para não poluir o gráfico (ex: 8 mais frequentes)

top_species <- dados %>%
  count(especie, sort = TRUE) %>%
  slice_max(order_by = n, n = 8) %>%
  pull(especie)


dados_subset <- dados %>%
  filter(especie %in% top_species)

# Plot
ggplot(dados_subset,
       aes(x = PC1_nutri, y = produt_z_g.ano, color = especie)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_y_continuous(trans = "log10") +   # << log to “decompress”
  facet_wrap(~ especie, scales = "free") + # free x and y scales per panel
  labs(title = "Species responses to the Nutrient gradient (PC1)",
       x = "Nutrient gradient (PC1_nutri)",
       y = "Biomass (z-score, log10 scale)") +
  theme_minimal()

# Set working directory
setwd("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/graficos")

# Save with ggsave
ggsave("species_response_nutrients.jpeg",
       width = 2000, height = 1600, units = "px", dpi = 300)

##### forest plot of species-specific slopes (productivity vs. PC1_clima #####

# ==== PACKAGES ====
library(readxl)
library(dplyr)
library(purrr)
library(broom)
library(ggplot2)
library(forcats)
library(stringr)

# ==== INPUTS ====
# Adjust paths and sheet name as needed:
excel_path <- "C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/dados_desempenho.xlsx"
sheet_name <- "principal"
out_dir    <- "C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/graficos"
min_n      <- 10   # increase to 20 for more stable CIs

# Column name for productivity (use your exact column name)
prod_col <- "produt_z_g.ano"

# ==== LOAD DATA ====
dat <- read_xlsx(excel_path, sheet = sheet_name) |>
  mutate(
    especie   = as.character(especie),
    PC1_nutri = as.numeric(PC1_nutri)
  )

stopifnot(prod_col %in% names(dat))

# ==== HELPER: fit slopes per species (simple lm: y ~ x) ====
fit_slopes <- function(df, response, predictor = "PC1_nutri", min_n = 10) {
  df |>
    select(especie, !!sym(predictor), !!sym(response)) |>
    rename(x = !!sym(predictor), y = !!sym(response)) |>
    filter(is.finite(x), is.finite(y)) |>
    group_by(especie) |>
    filter(n() >= min_n, var(x) > 0) |>      # ensure enough data & variation in x
    group_modify(~{
      m <- lm(y ~ x, data = .x)
      tb <- tidy(m)
      tibble(
        n       = nrow(.x),
        slope   = coef(m)[["x"]],
        lwr     = confint(m)["x", 1],
        upr     = confint(m)["x", 2],
        p_value = tb$p.value[tb$term == "x"]
      )
    }) |>
    ungroup() |>
    arrange(slope) |>
    mutate(especie = fct_inorder(especie))
}

# ==== HELPER: forest plot ====
plot_forest <- function(slopes_tbl, title, xlab) {
  ggplot(slopes_tbl,
         aes(x = slope, y = fct_reorder(especie, slope))) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2) +
    geom_point(size = 2) +
    labs(title = title,
         x = xlab,
         y = "Species") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white"),
      plot.background  = element_rect(fill = "white"),
      axis.text        = element_text(color = "black"),
      axis.title       = element_text(color = "black"),
      title            = element_text(color = "black")
    )
}

# ==== PRODUCTIVITY ====
slopes_prod <- fit_slopes(dat, response = prod_col, min_n = min_n)

forest_prod <- plot_forest(
  slopes_prod,
  title = "Species-specific responses to nutrients (PC1): Productivity",
  xlab  = "Slope of productivity (z-score) vs PC1_clima (β)"
)

# ==== SAVE OUTPUTS ====
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(out_dir, "forest_slopes_PC1nutri_productivity.jpeg"),
       forest_prod,
       width = 8,
       height = max(6, 0.4 * nrow(slopes_prod)),
       units = "in",
       dpi = 300)

write.csv(slopes_prod,
          file.path(out_dir, "forest_slopes_PC1nutri_productivity.csv"),
          row.names = FALSE)

