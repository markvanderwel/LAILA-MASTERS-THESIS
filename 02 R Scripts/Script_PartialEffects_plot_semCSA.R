############################
### Plot Partial Effects ###
############################

# Pacotes
library(lme4)
library(lmerTest)    # p-values para os efeitos fixos
library(ggeffects)   # efeitos marginais (ggpredict)
library(ggplot2)
library(patchwork)   # montar os painéis

# --- Model -----------------------------------------------------------
m12 <- lmer(log_produt ~ c.n_soloid + SR + pcps1ab + (1 | site),
            data = dadosmisto1, REML = FALSE)

# p-values dos efeitos fixos (para mostrar no subtítulo)
pv <- summary(m12)$coefficients[, "Pr(>|t|)"]
p_txt <- sprintf("p(c.n_soloid)=%.3f | p(SR)=%.3f | p(pcps1ab)=%.3f",
                 pv["c.n_soloid"], pv["SR"], pv["pcps1ab"])

# Efeitos marginais no nível populacional (sem efeitos aleatórios)
eff_c   <- ggpredict(m12, terms = "c.n_soloid", type = "fixed")
eff_sr  <- ggpredict(m12, terms = "SR",         type = "fixed")
eff_pc1 <- ggpredict(m12, terms = "pcps1ab",    type = "fixed")

# --- Função para plotar pontos + linha predita (escala log) -------------------
plot_log_pts <- function(eff, xvar, xlab){
  ggplot() +
    # pontos observados
    geom_point(data = dadosmisto1, aes_string(x = xvar, y = "log_produt"),
               color = "gray40", size = 1.8, alpha = 0.6) +
    # faixa de confiança do modelo
    geom_ribbon(data = eff, aes(x = x, ymin = conf.low, ymax = conf.high),
                fill = "skyblue3", alpha = 0.25) +
    # linha predita
    geom_line(data = eff, aes(x = x, y = predicted),
              color = "blue4", linewidth = 1.1) +
    labs(x = xlab, y = "Predito: log(produt)") +
    theme_minimal(base_size = 12)
}

# --- Função para plotar pontos + linha predita (escala original) --------------
plot_orig_pts <- function(eff, xvar, xlab){
  ggplot() +
    geom_point(data = dadosmisto1,
               aes_string(x = xvar, y = "exp(log_produt)"),
               color = "gray40", size = 1.8, alpha = 0.6) +
    geom_ribbon(data = eff,
                aes(x = x, ymin = exp(conf.low), ymax = exp(conf.high)),
                fill = "skyblue3", alpha = 0.25) +
    geom_line(data = eff, aes(x = x, y = exp(predicted)),
              color = "blue4", linewidth = 1.1) +
    labs(x = xlab, y = "Predito: produt (escala original)") +
    theme_minimal(base_size = 12)
}

# --- Montar painéis -----------------------------------------------------------
g1_log <- plot_log_pts(eff_c,   "c.n_soloid", "C:N do solo (c.n_soloid)")
g2_log <- plot_log_pts(eff_sr,  "SR",         "Riqueza de espécies (SR)")
g3_log <- plot_log_pts(eff_pc1, "pcps1ab",    "PCPS1 (pcps1ab)")

fig_log <- (g1_log | g2_log | g3_log) +
  plot_annotation(title = "Efeitos marginais com dados observados — escala log",
                  subtitle = p_txt)

# --- Também pode fazer na escala original -------------------------------------
g1_or <- plot_orig_pts(eff_c,   "c.n_soloid", "Soil C:N ratio")
g2_or <- plot_orig_pts(eff_sr,  "SR",         "Species Richness (SR)")
g3_or <- plot_orig_pts(eff_pc1, "pcps1ab",    "First PCPS axis")

fig_or <- (g1_or | g2_or | g3_or) +
  plot_annotation(title = "Partial effects with observed points — original scale",
                  subtitle = p_txt)

fig_or

setwd("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/graficos")
ggsave("fig_m12_efeitos_original_pts.png", fig_or, width = 11, height = 4.2, dpi = 300)


##########################################
##########################################
##########################################


# Pacotes
library(lme4)
library(lmerTest)
library(ggeffects)
library(dplyr)
library(ggplot2)
library(patchwork)

# --------------------------------------------------------------------
# Modelo (usa nomes originais do seu data frame)
# --------------------------------------------------------------------
m12 <- lmer(log_produt ~ c.n_soloid + SR + pcps1ab + (1 | site),
            data = dadosmisto1, REML = FALSE)

# p-values (para subtítulo)
pv <- summary(m12)$coefficients[, "Pr(>|t|)"]
p_txt_en <- sprintf("p(C:N)=%.3f | p(SR)=%.3f | p(PCPS1)=%.3f",
                    pv["c.n_soloid"], pv["SR"], pv["pcps1ab"])

# --------------------------------------------------------------------
# Data frame para PLOT com nomes em inglês (modelo continua com originais)
# --------------------------------------------------------------------
dados_plot <- dadosmisto1 %>%
  rename(
    soil_C_N_ratio      = c.n_soloid,
    species_richness_SR = SR,
    PCPS1               = pcps1ab,
    log_productivity    = log_produt
  )

# --------------------------------------------------------------------
# Efeitos marginais (somente efeitos fixos)
# --------------------------------------------------------------------
eff_c   <- ggpredict(m12, terms = "c.n_soloid", type = "fixed")
eff_sr  <- ggpredict(m12, terms = "SR",         type = "fixed")
eff_pc1 <- ggpredict(m12, terms = "pcps1ab",    type = "fixed")

# --------------------------------------------------------------------
# Função de plot (escala original) com pontos + linha + CI
# --------------------------------------------------------------------
base_sz <- 11
thin_margin <- margin(4, 6, 4, 6)

plot_orig_pts <- function(eff, df, xvar, xlab){
  ggplot() +
    geom_point(data = df,
               aes_string(x = xvar, y = "exp(log_productivity)"),
               color = "gray40", size = 1.9, alpha = 0.65) +
    geom_ribbon(data = eff,
                aes(x = x, ymin = exp(conf.low), ymax = exp(conf.high)),
                alpha = 0.18, fill = "#6baed6") +
    geom_line(data = eff, aes(x = x, y = exp(predicted)),
              linewidth = 1.15, color = "#08519c") +
    labs(x = xlab, y = "Predicted: productivity (original scale)") +
    theme_minimal(base_size = base_sz) +
    theme(plot.margin = thin_margin)
}

# Gráficos individuais
g_cn   <- plot_orig_pts(eff_c,   dados_plot, "soil_C_N_ratio",      "Soil C:N ratio")
g_sr   <- plot_orig_pts(eff_sr,  dados_plot, "species_richness_SR", "Species Richness (SR)") +
  ylab(NULL)  # remove o y-label do painel da direita
g_pcps <- plot_orig_pts(eff_pc1, dados_plot, "PCPS1",               "First PCPS axis")

library(patchwork)

# design: A (top-left), B (top-right); C (bottom-left), vazio à direita
design <- "
AB
C.
"

# MAPEIE PELOS NOMES!  👉 A=g_cn, B=g_sr, C=g_pcps
fig_or <- wrap_plots(
  A = g_cn,        # (a) Soil C:N  — topo esquerdo
  B = g_sr,        # (b) SR        — topo direito
  C = g_pcps,      # (c) PCPS1     — embaixo esquerdo
  design = design,
  heights = c(1, 0.55)   # deixa a linha de baixo mais baixa
) + plot_annotation(
  title = "Partial effects with observed points — original scale",
  subtitle = p_txt_en,
  tag_levels = "a",        # gera (a), (b), (c) exatamente nessa ordem A,B,C
  tag_prefix = "(", tag_suffix = ")"
)

fig_or
ggsave("fig_m12_partial_effects_original_EN.png", fig_or, width = 11, height = 6.4, dpi = 300)
