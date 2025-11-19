############################################
#### MASTER SCRIPT – DIVERSITY METRICS  ####
############################################

# 0) Load packages
library(dplyr)
library(tibble)
library(picante)
library(vegan)
library(readr)

# -------------------------------------------------
# 1) Input objects already in the environment
# -------------------------------------------------
# comunidade            -> community matrix (sites x species) for phylogenetic metrics
# pdis.ord              -> phylogenetic distance matrix (matching 'comunidade')
# tree_ok               -> phylogenetic tree (class 'phylo') for PD and psd()
# comunidade_diversity  -> community matrix used for Shannon / Simpson (can be the same as 'comunidade')

# Make sure sites are in the rownames:
# rownames(comunidade)
# rownames(comunidade_diversity)

# -------------------------
# 2) Faith's PD (picante)
# -------------------------
pd_faith <- pd(comunidade, tree_ok, include.root = FALSE)

pd_only <- pd_faith |>
  rownames_to_column("site") |>
  select(site, PD = PD)   # rename PD column if you want a clearer name

# -------------------------
# 3) MPD and MNTD (abundance-weighted)
# -------------------------
mpd_df <- tibble(
  site = rownames(comunidade),
  mpd  = as.numeric(mpd(comunidade, pdis.ord, abundance.weighted = TRUE))
)

mntd_df <- tibble(
  site = rownames(comunidade),
  mntd = as.numeric(mntd(comunidade, pdis.ord, abundance.weighted = TRUE))
)

# -----------------------------------------
# 4) SESPD, SESMPD, SESMNTD (z-scores)
# -----------------------------------------
ses_pd <- ses.pd(
  samp       = comunidade,
  tree       = tree_ok,
  null.model = "richness",
  runs       = 99
)

sespd_only <- ses_pd |>
  as.data.frame() |>
  rownames_to_column("site") |>
  select(site, sesPD = pd.obs.z)

ses_mpd <- ses.mpd(
  samp        = comunidade,
  dis         = pdis.ord,
  null.model  = "richness",
  runs        = 99
)

sesmpd_only <- ses_mpd |>
  as.data.frame() |>
  rownames_to_column("site") |>
  select(site, sesMPD = mpd.obs.z)

ses_mntd <- ses.mntd(
  samp        = comunidade,
  dis         = pdis.ord,
  null.model  = "richness",
  runs        = 99
)

sesmntd_only <- ses_mntd |>
  as.data.frame() |>
  rownames_to_column("site") |>
  select(site, sesMNTD = mntd.obs.z)

# -----------------------------------------
# 5) PSV, PSC, PSR, PSE (from psd)
# -----------------------------------------
psd_all <- psd(comunidade, tree_ok)

psd_df <- psd_all |>
  as.data.frame() |>
  rownames_to_column("site") |>
  select(site, PSV, PSC, PSE, PSR)

# -----------------------------------------
# 6) Shannon and unbiased Simpson
# -----------------------------------------
shannon_df <- tibble(
  site    = rownames(comunidade_diversity),
  shannon = as.numeric(
    diversity(comunidade_diversity, index = "shannon")
  )
)

unbias_simp_df <- tibble(
  site        = rownames(comunidade_diversity),
  unbias_simp = as.numeric(
    simpson.unb(comunidade_diversity)
  )
)

# -----------------------------------------
# 7) Join all metrics into a single data frame
# -----------------------------------------
data_diversity <- pd_only %>%
  left_join(mpd_df,         by = "site") %>%
  left_join(mntd_df,        by = "site") %>%
  left_join(sespd_only,     by = "site") %>%
  left_join(sesmpd_only,    by = "site") %>%
  left_join(sesmntd_only,   by = "site") %>%
  left_join(psd_df,         by = "site") %>%
  left_join(unbias_simp_df, by = "site") %>%
  left_join(shannon_df,     by = "site")

# Quick check
glimpse(data_diversity)

# -----------------------------------------
# 8) Save as CSV
# -----------------------------------------
write_csv(
  data_diversity,
  "~/01 Masters_LA/MASTERS-THESIS/02 Datasets/02 processed_data/data_diversity.csv"
)
