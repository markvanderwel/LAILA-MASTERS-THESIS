############################################################
### SPECIES-LEVEL PERFORMANCE ANALYSIS
### Productivity responses to environmental gradients
### Author: Laíla Arnauth
### Repository: Master’s Thesis – Functional & Phylogenetic Diversity
############################################################

# ============================
# 1. Load required packages
# ============================
library(readxl)
library(tidyverse)
library(writexl)
library(lmerTest)
library(broom)
library(forcats)
library(purrr)
library(ggplot2)

# ============================
# 2. Load dataset
# ============================
data_path <- "01 Datasets/01_raw_data/dados_desempenho_semcsa.xlsx"
dados <- read_excel(data_path)

# Mixed model example (not used below but kept for reference)
lmer(produt_z_g.ano ~ especie * PC1_nutri + (1 | site.parcela), data = dados)

# ============================
# 3. Select the most frequent species (top 8)
# ============================
top_species <- dados %>%
  count(especie, sort = TRUE) %>%
  slice_max(order_by = n, n = 8) %>%
  pull(especie)

dados_subset <- dados %>%
  filter(especie %in% top_species)

# ============================
# 4. Scatterplots of species responses to PC1
# ============================
p_response <- ggplot(
  dados_subset,
  aes(x = PC1_nutri, y = produt_z_g.ano, color = especie)
) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_y_continuous(trans = "log10") +  # decompress variation
  facet_wrap(~ especie, scales = "free") +
  labs(
    title = "Species responses to the nutrient gradient (PC1)",
    x = "Nutrient gradient (PC1_nutri)",
    y = "Biomass (z-score, log10 scale)"
  ) +
  theme_minimal()

# Save plot
ggsave(
  "species_response.jpeg",
  p_response,
  width = 2000,
  height = 1600,
  units = "px",
  dpi = 300
)

############################################################
### FOREST PLOT — species-specific slopes (productivity vs PC1_clima)
############################################################

# ============================
# 5. Inputs
# ============================
excel_path <- "01 Datasets/01_raw_data/dados_desempenho_semcsa.xlsx"
sheet_name <- "principal"

output_dir <- "~/01 Masters_LA/06 Figures/02 plots"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

min_n <- 10  # minimum number of points per species for slope estimation
prod_col <- "produt_z_g.ano"  # response variable

# ============================
# 6. Load main dataset
# ============================

dat <- read_xlsx(excel_path, sheet = sheet_name) %>%
  mutate(
    especie   = as.character(especie),
    PC1_clima = as.numeric(PC1_clima)
  )

stopifnot(prod_col %in% names(dat))

# ============================
# 7. Helper function:
#    Fit per-species slope (lm: y ~ x)
# ============================
fit_slopes <- function(df, response, predictor = "PC1_clima", min_n = 10) {
  
  # Extract only required columns (base R to avoid masking issues)
  df2 <- df[, c("especie", predictor, response)]
  names(df2) <- c("especie", "x", "y")
  
  # Remove non-finite values
  df2 <- df2[is.finite(df2$x) & is.finite(df2$y), ]
  
  # Split dataset by species
  by_species <- split(df2, df2$especie)
  
  # Fit linear models
  results <- lapply(by_species, function(d) {
    if (nrow(d) < min_n) return(NULL)
    if (var(d$x) == 0)    return(NULL)  # cannot compute slope
    
    m  <- lm(y ~ x, data = d)
    ci <- suppressMessages(confint(m))
    tb <- suppressMessages(broom::tidy(m))
    
    data.frame(
      especie = unique(d$especie),
      n       = nrow(d),
      slope   = coef(m)[["x"]],
      lwr     = ci["x", 1],
      upr     = ci["x", 2],
      p_value = tb$p.value[tb$term == "x"],
      stringsAsFactors = FALSE
    )
  })
  
  out <- do.call(rbind, results)
  
  # Order species by slope
  out <- out[order(out$slope), ]
  out$especie <- factor(out$especie, levels = out$especie)
  
  return(out)
}

# ============================
# 8. Fit slopes for productivity
# ============================
slopes_prod <- fit_slopes(dat, response = prod_col, min_n = min_n)

# ============================
# 9. Forest plot function
# ============================
plot_forest <- function(slopes_tbl, title, xlab) {
  ggplot(slopes_tbl,
         aes(x = slope, y = fct_reorder(especie, slope))) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2) +
    geom_point(size = 2) +
    labs(title = title, x = xlab, y = "Species") +
    theme_minimal(base_size = 12)
}

forest_prod <- plot_forest(
  slopes_prod,
  title = "Species-specific responses to aridity (PC1): Productivity",
  xlab  = "Slope of productivity (z-score) vs PC1_clima (β)"
)

# ============================
# 10. Save outputs
# ============================
ggsave(
  file.path(output_dir, "forest_slopes_PC1_productivity.jpeg"),
  forest_prod,
  width = 8,
  height = max(6, 0.4 * nrow(slopes_prod)),
  units = "in",
  dpi = 300
)

# ============================
# 11. Keep only significant species (p < 0.05)
# ============================
slopes_sig <- slopes_prod %>%
  filter(!is.na(p_value), p_value < 0.05)

# If you want to check how many:
nrow(slopes_sig)

# Forest plot with significant slopes only
forest_sig <- plot_forest(
  slopes_sig,
  title = "Species with significant responses to aridity (PC1): Productivity",
  xlab  = "Slope of productivity (z-score) vs PC1_clima (β)"
)

# Save plot
ggsave(
  file.path(output_dir, "forest_slopes_PC1_productivity_significant.jpeg"),
  forest_sig,
  width = 8,
  height = max(6, 0.4 * nrow(slopes_sig)),
  units = "in",
  dpi = 300
)


############################################################
### SPECIES-LEVEL PERFORMANCE — NUTRIENT GRADIENT (PC1_nutri)
### Author: Laíla Arnauth
### Repository: Master’s Thesis – Functional & Phylogenetic Diversity
############################################################

# ============================
# 1. Load required packages
# ============================
library(readxl)
library(tidyverse)
library(writexl)
library(lmerTest)
library(broom)
library(forcats)
library(purrr)
library(ggplot2)

# ============================
# 2. Paths and general settings
# ============================
# All paths are relative to the project root (recommended for GitHub)
data_scatter_path <- "01 Datasets/01_raw_data/dados_desempenho_semcsa.xlsx"
data_forest_path  <- "01 Datasets/01_raw_data/dados_desempenho.xlsx"
output_dir        <- "06 Figures/02 plots"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

############################################################
### PART A — Scatterplots of species responses to PC1_nutri
############################################################

# ============================
# 3. Load data for scatterplots
# ============================
dados <- read_excel(data_scatter_path)

# Example mixed model (kept for reference)
lmer(
  produt_z_g.ano ~ especie * PC1_nutri + (1 | site.parcela),
  data = dados
)

# ============================
# 4. Select most frequent species (top 8)
# ============================
top_species <- dados %>%
  count(especie, sort = TRUE) %>%
  slice_max(order_by = n, n = 8) %>%
  pull(especie)

dados_subset <- dados %>%
  filter(especie %in% top_species)

# ============================
# 5. Scatterplot: productivity vs nutrient gradient
# ============================
p_nutrients <- ggplot(
  dados_subset,
  aes(x = PC1_nutri, y = produt_z_g.ano, color = especie)
) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_y_continuous(trans = "log10") +  # decompress variation
  facet_wrap(~ especie, scales = "free") +
  labs(
    title = "Species responses to the nutrient gradient (PC1_nutri)",
    x     = "Nutrient gradient (PC1_nutri)",
    y     = "Biomass (z-score, log10 scale)"
  ) +
  theme_minimal()

# Save scatterplot
ggsave(
  filename = file.path(output_dir, "species_response_nutrients.jpeg"),
  plot     = p_nutrients,
  width    = 2000,
  height   = 1600,
  units    = "px",
  dpi      = 300
)

############################################################
### PART B — Forest plot of species-specific slopes
###          Productivity vs PC1_nutri
############################################################

# ============================
# 6. Inputs for forest plot
# ============================
sheet_name <- "principal"
min_n      <- 10                     # minimum n per species
prod_col   <- "produt_z_g.ano"       # response variable

# ============================
# 7. Load data
# ============================
dat <- read_xlsx(data_forest_path, sheet = sheet_name) %>%
  mutate(
    especie   = as.character(especie),
    PC1_nutri = as.numeric(PC1_nutri)
  )

stopifnot(prod_col %in% names(dat))

# ============================
# 8. Helper: fit slopes per species (lm: y ~ x)
#    (base R inside to avoid select()/rename() conflicts)
# ============================
fit_slopes <- function(df, response, predictor = "PC1_nutri", min_n = 10) {
  
  # Extract only required columns
  df2 <- df[, c("especie", predictor, response)]
  names(df2) <- c("especie", "x", "y")
  
  # Remove non-finite values
  df2 <- df2[is.finite(df2$x) & is.finite(df2$y), ]
  
  # Split by species
  by_species <- split(df2, df2$especie)
  
  # Fit one model per species
  results <- lapply(by_species, function(d) {
    if (nrow(d) < min_n) return(NULL)
    if (var(d$x) == 0)    return(NULL)  # no variation in predictor
    
    m  <- lm(y ~ x, data = d)
    ci <- suppressMessages(confint(m))
    tb <- suppressMessages(broom::tidy(m))
    
    data.frame(
      especie = unique(d$especie),
      n       = nrow(d),
      slope   = coef(m)[["x"]],
      lwr     = ci["x", 1],
      upr     = ci["x", 2],
      p_value = tb$p.value[tb$term == "x"],
      stringsAsFactors = FALSE
    )
  })
  
  out <- do.call(rbind, results)
  
  # Order species by slope
  out <- out[order(out$slope), ]
  out$especie <- factor(out$especie, levels = out$especie)
  
  return(out)
}

# ============================
# 9. Helper: forest plot
# ============================
plot_forest <- function(slopes_tbl, title, xlab) {
  ggplot(
    slopes_tbl,
    aes(x = slope, y = fct_reorder(especie, slope))
  ) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2) +
    geom_point(size = 2) +
    labs(
      title = title,
      x     = xlab,
      y     = "Species"
    ) +
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

# ============================
# 10. Fit slopes for productivity
# ============================
slopes_prod <- fit_slopes(dat, response = prod_col, min_n = min_n)

# Forest plot with all species
forest_prod <- plot_forest(
  slopes_prod,
  title = "Species-specific responses to nutrient gradient (PC1): Productivity",
  xlab  = "Slope of productivity (z-score) vs PC1_nutri (β)"
)

# Save full forest plot
ggsave(
  file.path(output_dir, "forest_slopes_PC1nutri_productivity.jpeg"),
  forest_prod,
  width  = 8,
  height = max(6, 0.4 * nrow(slopes_prod)),
  units  = "in",
  dpi    = 300
)


# ============================
# 11. Significant species only (p < 0.05)
# ============================
slopes_sig <- slopes_prod %>%
  filter(!is.na(p_value), p_value < 0.05)

forest_sig <- plot_forest(
  slopes_sig,
  title = "Species with significant responses to nutrient gradient (PC1): Productivity",
  xlab  = "Slope of productivity (z-score) vs PC1_nutri (β)"
)

# Save significant-only forest plot
ggsave(
  file.path(output_dir, "forest_slopes_PC1nutri_productivity_significant.jpeg"),
  forest_sig,
  width  = 8,
  height = max(6, 0.4 * nrow(slopes_sig)),
  units  = "in",
  dpi    = 300
)


