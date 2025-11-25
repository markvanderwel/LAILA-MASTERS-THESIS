
########################################
### SCRIPT FOR PHYLOGENETIC ANALYSIS ###
########################################

# Last update: 2025/11/19 (YYYY/MM/DD)
# Authors: Laíla Arnauth & André T. C. Dias
# Post Graduate Program in Ecology - UFRJ - Brazil
# Collaboration with University of Regina - Canada # Forest Dynamics Lab


# Load all packages

library(readr)
library(tidyverse)
library(readxl)
library(writexl)
library(ggplot2)
library(corrplot)
library("V.PhyloMaker2")
library(picante)
library(ape)
library(vegan)
library(lme4)
library(lmerTest) # LMM
library(MuMIn) #AICc
library(car)
library(broom.mixed)
library(grid)  # unit()
library(sjPlot)
library(ggeffects)   # ggpredict
library(patchwork)
library(effects)
library(PCPS)
library(FactoMineR)  # advanced PCA
library(factoextra)
library(hillR) # Hill numbers
library(phytools)

#########################
### PHYLOGENETIC TREE ###
#########################

# Install V.PhyloMaker2
# install.packages("devtools")
# devtools::install_github("jinyizju/V.PhyloMaker2")

# input the sample species list
example <- read_csv("~/01 Masters_LA/00 MASTERS-DATA/01 Datasets/01_raw_data/filogenia_total.csv")
View(example)

### generate a phylogeny for the sample species list
tree <- phylo.maker(example, tree = GBOTB.extended.TPL,output.sp.list = TRUE,nodes = nodes.info.1.TPL, scenarios="S3")
tree_ok <- tree$scenario.3
summary(tree_ok)

# Plot Phylo Tree

plot.phylo(tree_ok, type = "tidy", show.tip.label = TRUE, show.node.label = FALSE,edge.color = "gray",edge.width = 1,tip.color = "black")

png("~/01 Masters_LA/06 Figures/03 Plot_phylo_trees/arvore_filog.jpeg",width = 12, height = 10, units = "in", res = 300)
plot.phylo(tree_ok, 
           type = "tidy", 
           show.tip.label = TRUE, 
           show.node.label = FALSE, 
           edge.color = "gray", 
           edge.width = 2, 
           tip.color = "black", 
           cex = 0.7,
           direction = "upward") 
dev.off()

# The 'cophenetic' function calculates dissimilarity based on the phylogenetic distances of the tree
phylo.dissim = cophenetic(tree_ok)
View(phylo.dissim)

# We now standardize distances on a scale ranging from 0 to 1.
std.pdis = (phylo.dissim-min(phylo.dissim))/(max(phylo.dissim)-min(phylo.dissim))
View(std.pdis)

# We arranged the species in alphabetical order, as in the community matrix.
pdis.ord = std.pdis[order(row.names(std.pdis)), order(row.names(std.pdis))]
View(pdis.ord)

### Importing the community matrix
## comunidade_pd uses Basal Area, and not abundance; basal area is used in the calculation of Biomass, so it's prefereable to use Abundance

comunidade <- read.csv("01 Datasets/01_raw_data/comunidade_abund.csv",
  row.names = 1,
  header = TRUE,
  sep = ";")

### ------------------------------------------------ ###
### phylogenetic diversity is calculated w/ picante  ###
### ------------------------------------------------ ###

colnames(comunidade)==colnames(pdis.ord) # checking if they match

# PD FAITH

pd.faith <- pd(comunidade,tree_ok,include.root = FALSE)
view(pd.faith)
write_xlsx(pd.faith, "~/01 Masters_LA/00 MASTERS-DATA/01 Datasets/02_processed_data/PD_Faith.xlsx")

pd_only <- pd.faith %>%
  tibble::rownames_to_column(var = "site") %>%
  dplyr::select(-SR)


# MPD

mpd_df <- tibble::tibble(
  site = rownames(comunidade),
  mpd  = as.numeric(mpd(comunidade, pdis.ord, abundance.weighted = TRUE))
)
write_xlsx(mpd_df, "~/01 Masters_LA/00 MASTERS-DATA/01 Datasets/02_processed_data/MPD.xlsx")

# MNTD (species are more broadly or finely spaced in a community compared to another?)

mntd_df <- tibble::tibble(
  site = rownames(comunidade),
  mntd = as.numeric(mntd(comunidade, pdis.ord, abundance.weighted = TRUE))
)

write_xlsx(mntd_df, "~/01 Masters_LA/00 MASTERS-DATA/01 Datasets/02_processed_data/MNTD.xlsx")

### ---------------------------------- ###
### Standardized effect size with ape  ###
### ---------------------------------- ###

### Standardized effect size of PD  ##### Runs = 999

# tree_clean <- drop.tip(tree_ok, tip = which(is.na(tree_ok$tip.label))) #OPTIONAL

ses_pd <- ses.pd(comunidade, tree_ok, null.model="richness", runs = 999,iterations = 1000,include.root = FALSE)
ses_pd <- as.data.frame(ses_pd)

write_xlsx(ses_pd, "~/01 Masters_LA/00 MASTERS-DATA/01 Datasets/02_processed_data/ses_pd.xlsx")

ses_pd <- ses_pd |> 
  tibble::rownames_to_column(var = "site")

sespd_only <- ses_pd |> 
  transmute(site, sespd = pd.obs.z)
view(sespd_only)

### Standardized effect size of MPD ##### Runs = 999

ses_mpd <- ses.mpd(comunidade, pdis.ord,null.model="richness")
ses_mpd <- as.data.frame(ses_mpd)

write_xlsx(ses_mpd, "~/01 Masters_LA/00 MASTERS-DATA/01 Datasets/02_processed_data/ses_mpd.xlsx")

ses_mpd <- ses_mpd |> 
  tibble::rownames_to_column(var = "site")

sesmpd_only <- ses_mpd |> 
  transmute(site, sesmpd = mpd.obs.z)
view(sesmpd_only)

### Standardized effect size of MNTD ##### Runs = 999

ses_mntd <- ses.mntd(comunidade, pdis.ord,null.model="richness")
ses_mntd <- as.data.frame(ses_mntd)
write_xlsx(ses_mntd, "~/01 Masters_LA/00 MASTERS-DATA/01 Datasets/02_processed_data/ses_mntd.xlsx")

ses_mntd <- ses_mntd |> 
  tibble::rownames_to_column(var = "site")

sesmntd_only <- ses_mntd |> 
  transmute(site, sesmntd = mntd.obs.z)
view(sesmntd_only)

#######
# PSV #
####### picante

#'psv(samp,tree,compute.var=TRUE,scale.vcv=TRUE)
#'psr(samp,tree,compute.var=TRUE,scale.vcv=TRUE)
#'pse(samp,tree,scale.vcv=TRUE)
#'psc(samp,tree,scale.vcv=TRUE)
#'psd(samp,tree,compute.var=TRUE,scale.vcv=TRUE)
#'psv.spp(samp,tree)

# Phylogenetic Species variability: quantifies how phylogenetic relatedness decreases the variance of a hypothetical unselected/neutral trait shared by all species in a community.

# PSD function: all at once

PSD <- psd(comunidade,tree_ok)
view(PSD)
write_xlsx(PSD, "~/01 Masters_LA/MASTERS-THESIS/02 Datasets/02 processed_data/PSD.xlsx")

psd_df <- PSD %>%
  rownames_to_column("site") %>%
  select(site, PSV, PSC, PSE, PSR)
view(psd_df)

### OPTIONAL
# Phylogenetic Species Richness: PSV * SR
PSR <- psr(comunidade,tree_ok)
write_xlsx(PSR, "~/01 Masters_LA/00 MASTERS-DATA/01 Datasets/02_processed_data/PSR.xlsx")

# Phylogenetic Species Evenness: is the metric PSV modified to incorporate relative species abundances
PSE <- pse(comunidade,tree_ok)
write_xlsx(PSE, "~/01 Masters_LA/00 MASTERS-DATA/01 Datasets/02_processed_data/PSE.xlsx")

# Phylogenetic species clustering: is a metric of the branch tip clustering of species across a sample’s phylogeny

PSC <- psc(comunidade,tree_ok)
write_xlsx(PSC, "~/01 Masters_LA/00 MASTERS-DATA/01 Datasets/02_processed_data/PSC.xlsx")

###########################
### TAXONOMIC DIVERSITY ###
###########################

comunidade_diversity <- read.csv("~/01 Masters_LA/00 MASTERS-DATA/01 Datasets/01_raw_data/comunidade_diversity.csv",
                                 row.names = 1,
                                 header = TRUE,
                                 sep = ";")
# Simpson
unbias_simp_df <- tibble::tibble(
  site = rownames(comunidade_diversity),
  unbias_simp = as.numeric(simpson.unb(comunidade_diversity))
)
write_xlsx(unbias.simp.df, "~/01 Masters_LA/00 MASTERS-DATA/01 Datasets/02_processed_data/Unbias_Simp.xlsx")

#Shannon
shannon_df <- tibble::tibble(
  site = rownames(comunidade_diversity),
  shannon = as.numeric(diversity(comunidade_diversity, index = "shannon"))
)

write_xlsx(shannon_df, "~/01 Masters_LA/00 MASTERS-DATA/01 Datasets/02_processed_data/shannon.xlsx")

### Saving diversity data with community weighted by abundance

data_diversity <- pd_only %>%
  left_join(mpd_df,        by = "site") %>%
  left_join(mntd_df,       by = "site") %>%
  left_join(sespd_only,    by = "site") %>%
  left_join(sesmpd_only,   by = "site") %>%
  left_join(sesmntd_only,  by = "site") %>%
  left_join(psd_df,        by = "site") %>%
  left_join(unbias_simp_df, by = "site") %>%
  left_join(shannon_df,    by = "site")

write_csv(
  data_diversity,
  "01 Datasets/02_processed_data/data_diversity.csv"
)

write_xlsx(data_diversity,
           "01 Datasets/02_processed_data/data_diversity.xlsx")

########################
### DATA EXPLORATION ###
########################

dadosmisto <- read.csv("01 Datasets/01_raw_data/dadosmisto.csv",
                                 header = TRUE)

dadosmisto1 <- dadosmisto[-c(1:7), ]

png("~/01 Masters_LA/06 Figures/01 exploratory_plots/hist_biomass.png", width = 800, height = 600)
hist(dadosmisto$biomassa_z_kg, main="Histograma Biomassa", xlab="Valores", col="blue", border="black")
dev.off()

png("~/01 Masters_LA/06 Figures/01 exploratory_plots/hist_logprodut.png", width = 800, height = 600)
hist(dadosmisto$log_produt, main="Histograma Log Produtividade", xlab="Valores", col="green", border="black")
dev.off()

# SHAPIRO TEST

shap_biomass <- shapiro.test(dadosmisto$biomassa_z_kg) # NOT NORMAL
dados_transformados <- log(dadosmisto$biomassa_z_kg + 1)  # Add 1 to avoid log(0)
dados_transformados.df <- as.data.frame(dados_transformados)
write_xlsx(dados_transformados.df, "01 Datasets/02_processed_data/logbiomassa.xlsx")
shapiro.test(dados_transformados) 

shap_produt <- shapiro.test(dadosmisto$produt_z_g.ano) # NÃO É NORMAL
dados_transformados <- log(dadosmisto$produt_z_g.ano + 1)  # Adicione 1 para evitar log(0)
dados_transformados.df <- as.data.frame(dados_transformados)
shapiro.test(dados_transformados) 
write_xlsx(dados_transformados.df, "01 Datasets/02_processed_data/logprodut.xlsx")


# DENSITY

plot(density(dadosmisto$log_biomass), main="Densidade dos Dados")
plot(density(dadosmisto$log_produt), main="Densidade dos Dados")

# CHECKING FOR OUTLIERS

png("~/01 Masters_LA/06 Figures/01 exploratory_plots/boxplot_LogBiomassa.png",
    width = 800, height = 600)
boxplot(dadosmisto$log_biomass, 
        main = "Boxplot - Log Biomass (kg)",
        col = "#8B008B")
dev.off()

png("~/01 Masters_LA/06 Figures/01 exploratory_plots/boxplot_LogProdutividade.png",
    width = 800, height = 600)
boxplot(dadosmisto$log_produt, 
        main = "Boxplot - Productivity (g/cm²/year)", 
        col = "#32CD32")
dev.off()

###################################
####### CORRELATION MATRIX ########
#### INDEPENDENT VARIABLES (X) ####
###################################

# Read dataset
data <- read.csv("01 Datasets/01_raw_data/dadosmisto.csv", header = TRUE)

# Remove first 7 rows (without CSA)
data_filtered <- data[-c(1:7), ]

# Define independent variables
labels_pretty <- c(
  ppt         = "Precipitation (ppt)",
  tmax        = "Max Temperature (°C)",
  tmin        = "Min Temperature (°C)",
  pet         = "PET",
  vpd         = "VPD",
  mcwd        = "MCWD",
  PC1_clima   = "Climate PC1",
  altitude    = "Altitude (m)",
  declividade = "Slope (%)",
  silte       = "Silt (%)",
  ph          = "Soil pH",
  valor_s     = "Soil S",
  valor_t     = "Soil T",
  PC1nutri    = "Nutrients PC1",
  PC2nutri    = "Nutrients PC2",
  season_temp = "Temp Seasonality",
  season_ppt  = "PPT Seasonality",
  c.n_soloid  = "Soil C/N Ratio"
)

variables_x <- c("ppt", "tmax", "tmin", "pet", "vpd", "mcwd", "PC1_clima",
                 "altitude", "declividade", "silte", "ph", "valor_s", "valor_t",
                 "PC1nutri", "PC2nutri", "season_temp", "season_ppt", "c.n_soloid")

# Create dataframe with X variables
data_x <- data_filtered %>%
  dplyr::select(all_of(variables_x))

# Correlation matrix
cor_matrix_x <- cor(data_x, use = "complete.obs")

# Replace row and column names with pretty labels
rownames(cor_matrix_x) <- labels_pretty[variables_x]
colnames(cor_matrix_x) <- labels_pretty[variables_x]

# Plot
png("~/01 Masters_LA/06 Figures/01 exploratory_plots/correlation_matrix_environmental_without_CSA.png",
    width = 1200, height = 800, res = 150)

corrplot(cor_matrix_x,
         method = "circle",
         type = "upper",
         tl.col = "black",
         tl.cex = 0.8)

dev.off()

##########################
### LINEAR MIXED MODEL ###
##########################

dadosmisto <- read.csv("01 Datasets/01_raw_data/dadosmisto.csv",
                       header = TRUE,
                       sep = ";")

dadosmisto1 <- dadosmisto[-c(1:7), ]

# OFICIAL ### site as random variable

modelo_nulo <- lmer(log_produt ~ (1|site), data = dadosmisto1, REML = FALSE)

### WITHOUT CSA ###

m1 <- lmer(log_produt ~ sr + silte + (1 | site), data = dadosmisto1, REML = FALSE)
summary(m1) 
anova(modelo_nulo, m1) # p-value = 0.005454 **
r.squaredGLMM(m1) # 0.5561193
AICc(m1) # 73.76227

m2 <- lmer(log_produt ~ sr + silte + pcps1 + (1 | site), data = dadosmisto1, REML = FALSE)
summary(m2) 
anova(modelo_nulo, m2) # p-value = 0.00144 **
r.squaredGLMM(m2) # 0.5609249
AICc(m2) # 71.12221

m3 <- lmer(log_produt ~ sr + pcps1 + (1 | site), data = dadosmisto1, REML = FALSE)
summary(m3) 
anova(modelo_nulo, m3) # p-value = 0.00109 **
r.squaredGLMM(m3) # 0.5643258
AICc(m3) # 70.54185

m4 <- lmer(log_produt ~ sespd + (1 | site), data = dadosmisto1, REML = FALSE)
summary(m4) 
anova(modelo_nulo, m4) # p-value = 0.9531    
r.squaredGLMM(m4) # 0.5117967
AICc(m4) # 81.83127

m5 <- lmer(log_produt ~ sr * pcps1 + (1 | site), 
           data = dadosmisto1, REML = FALSE)
summary(m5)
AICc(m5) # 72.37014

m6 <- lmer(log_produt ~ riq_inicial + (1 | site), data = dadosmisto1, REML = FALSE)
summary(m6) # p-value = 0.904   
r.squaredGLMM(m6) # 0.5127016
AICc(m6) # 81.81974

m7 <- lmer(log_produt ~ ph + (1 | site), data = dadosmisto1, REML = FALSE)
summary(m7) # p-value = 0.977
r.squaredGLMM(m7) # 0.51279
AICc(m7) # 81.83388

m8 <- lmer(log_produt ~ n_solo + (1 | site), data = dadosmisto1, REML = FALSE)
summary(m8) # p-value = 0.133     
r.squaredGLMM(m8) # 0.5010822
AICc(m8) # 79.57211

m9 <- lmer(log_produt ~ c.n_soloid + (1 | site), data = dadosmisto1, REML = FALSE)
summary(m9) # p-value =  0.0294 * 
r.squaredGLMM(m9) # 0.5271506
AICc(m9) # 76.89054

m10 <- lmer(log_produt ~ c.n_soloid + sr + ph + silte + pcps1 + (1 | site), data = dadosmisto1, REML = FALSE)
summary(m10) 
anova(modelo_nulo, m10) # p-value = 0.001474 **   
r.squaredGLMM(m10) # 0.5572333
AICc(m10) # 72.12193

m11 <- lmer(log_produt ~ c.n_soloid + sr + silte + pcps1 + (1 | site), data = dadosmisto1, REML = FALSE)
summary(m11) 
anova(modelo_nulo, m11) # p-value = 0.0005941 ***   
r.squaredGLMM(m11) # 0.5572206
AICc(m11) # 69.51542

m12 <- lmer(log_produt ~ c.n_soloid + sr + pcps1 + (1 | site), data = dadosmisto1, REML = FALSE)
summary(m12) 
anova(modelo_nulo, m12) # p-value = 0.0003766 *** 
r.squaredGLMM(m12) # 0.5667528
AICc(m12) # 68.28968
### m12 best model

m13 <- lmer(log_produt ~ c_soloid + (1 | site), data = dadosmisto1, REML = FALSE)
summary(m13) 
anova(modelo_nulo, m13) # p-value = 0.577
r.squaredGLMM(m13) # 0.5142692
AICc(m13) # 81.51027

m14 <- lmer(log_produt ~ season_ppt + (1 | site), data = dadosmisto1, REML = FALSE)
summary(m14) # p-value = 0.0387 * 
r.squaredGLMM(m14) # 0.5231383
AICc(m14) # 76.82315

m15 <- lmer(log_produt ~ season_temp + (1 | site), data = dadosmisto1, REML = FALSE)
summary(m15) # p-value = 0.408
r.squaredGLMM(m15) # 0.5237957
AICc(m15) # 80.57329

m16 <- lmer(log_produt ~ c.n_soloid + sr + pcps1 + season_ppt + (1 | site), data = dadosmisto1, REML = FALSE)
summary(m16) 
anova(modelo_nulo, m16) # p-value = error
r.squaredGLMM(m16) # 0.5667528
AICc(m16) # 66.55721

m17 <- lmer(log_produt ~ sr + pcps1 + season_ppt + (1 | site), data = dadosmisto1, REML = FALSE)
summary(m17) 
anova(modelo_nulo, m17) # p-value = error
r.squaredGLMM(m17) # 0.5741846
AICc(m17) # 66.78168
### NEW BEST MODEL

##### MULTICOLINEARITY ########

# install.packages("car")

m12 <- lmer(log_produt ~ c.n_soloid + sr + pcps1 + (1 | site), data = dadosmisto1, REML = FALSE)

# Calculando o VIF
vif(m12) # c.n_soloid   sr         pcps1 
         # 1.027845   1.117701   1.089277 

### CARBON ###

# Calcular conteúdo de carbono como 46.7% da produtividade (biomassa)
dadosmisto1$carbono <- dadosmisto1$produt_z_g.ano * 0.467

# Aplicar log no conteúdo de carbono
dadosmisto1$log_carbono <- log(dadosmisto1$carbono)

# MODELO NULO (sem preditores fixos)
modelo_nulo <- lmer(log_carbono ~ (1 | site), data = dadosmisto1, REML = FALSE)
summary(modelo_nulo)
AICc(modelo_nulo) # 86.88589

# MODELO COM VARIÁVEIS PREDITIVAS
c1 <- lmer(log_carbono ~ c.n_soloid + sr + silte + pcps1 + (1 | site), data = dadosmisto1, REML = FALSE)
summary(c1)
anova(modelo_nulo, c1)   # P-value = 0.0005941 ***
r.squaredGLMM(c1)        # R²c = 0.5572206
AICc(c1)                 # 69.51542

#############
### PLOTS ###
#############

## =============== ##
##  Coefplot m11   ##  
## =============== ##

dadosmisto <- read.csv("01 Datasets/01_raw_data/dadosmisto.csv",
                       header = TRUE,
                       sep = ";")

dadosmisto1 <- dadosmisto[-c(1:7), ]

## 1) Standardizing variables (mean = 0, sd = 1)

dados_std <- dadosmisto1
dados_std$c.n_soloid <- scale(dadosmisto1$c.n_soloid)
dados_std$SR         <- scale(dadosmisto1$sr)
dados_std$silte      <- scale(dadosmisto1$silte)
dados_std$pcps1ab    <- scale(dadosmisto1$pcps1)


## 2) Run standardized model

m12_std <- lmer(log_produt ~ c.n_soloid + SR + pcps1ab + (1 | site),
                data = dados_std, REML = FALSE)

summary(m12_std)

## 3) Extract coefficients + IC95% of the stand. model
coef_df <- tidy(m11_std, effects = "fixed", conf.int = TRUE)

## 4) Remove intercept and add significant p-value
coef_df_sem_intercepto <- subset(coef_df, term != "(Intercept)")
coef_df_sem_intercepto$significativo <- ifelse(coef_df_sem_intercepto$p.value < 0.05,
                                               "Significant", "Non-significant")

## 5) Labels and desirable order (SR → PCPS1 → C:N → Silt)
coef_df_sem_intercepto$term <- factor(coef_df_sem_intercepto$term,
                                      levels = c("SR","pcps1ab","c.n_soloid","silte"),
                                      labels = c("Species richness (SR)",
                                                 "Phylogenetic composition (PCPS1)",
                                                 "Soil C:N ratio",
                                                 "Silt (%)")
)
coef_df_sem_intercepto$term <- factor(coef_df_sem_intercepto$term,
                                      levels = rev(levels(coef_df_sem_intercepto$term))  # SR no topo
)

## 6) Coefplot (bolas pretas = significativo)
p <- ggplot(coef_df_sem_intercepto,
            aes(x = estimate, y = term, fill = significativo)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.22, color = "black") +
  geom_point(size = 3, shape = 21, color = "black", stroke = 1) +
  scale_fill_manual(values = c("Significant" = "black",
                               "Non-significant" = "white")) +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold"),
        legend.position = "none",
        plot.margin = unit(c(6, 10, 6, 10), "mm")) +
  labs(x = expression("Standardized regression coefficients ("*beta*" ± 95% CI)"),
       y = "Predictors")
p
## 7) Salvar sem cortar nada

ggsave("~/01 Masters_LA/06 Figures/02 plots/graf_coefplot_msel.jpeg",
       plot = p, width = 14, height = 9, units = "cm",
       dpi = 600, bg = "white", limitsize = FALSE)



################
### PLOT GLM ###
################

dadosmisto <- read.csv("01 Datasets/01_raw_data/dadosmisto.csv",
                       header = TRUE,
                       sep = ";")

dadosmisto1 <- dadosmisto[-c(1:7), ]


## SR 

fit <- lm(log_produt ~ sr, data = dadosmisto1)
summary(fit)

r2  <- summary(fit)$r.squared
p   <- coef(summary(fit))["sr", "Pr(>|t|)"]

subtxt <- sprintf("R\u00B2 = %.2f; p = %s",
                  r2,
                  ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))


g <- ggplot(dadosmisto1, aes(sr, log_produt)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, colour = "yellow") +
  labs(
    x = "Species richness (SR)",
    y = expression("log Biomass accumulation (g m"^-2*" year"^-1*")"),
    subtitle = subtxt
  ) +
  theme_minimal(base_size = 11)
g

ggsave("~/01 Masters_LA/06 Figures/02 plots/logprodut_SR.jpeg", g, width = 15, height = 10, units = "cm", dpi = 600)



## Silte

fit <- lm(log_produt ~ silte, data = dadosmisto1)
s   <- summary(fit)

r2  <- s$r.squared

coefs <- coef(s)
slope_row <- setdiff(rownames(coefs), "(Intercept)")
p  <- coefs[slope_row, "Pr(>|t|)"]

subtxt <- sprintf("R\u00B2 = %.2f; p %s",
                  r2,
                  ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

g <- ggplot(dadosmisto1, aes(silte, log_produt)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, colour = "brown") +
  labs(
    x = "Soil silt content (%)",
    y = expression("log Biomass accumulation (g m"^-2*" year"^-1*")"),
    subtitle = subtxt
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.margin = unit(c(6,10,6,10), "mm"))

g

ggsave("~/01 Masters_LA/06 Figures/02 plots/logprodut_silt.jpeg", g,
       width = 15, height = 10, units = "cm", dpi = 600, bg = "white")

## Soil C:N

fit <- lm(log_produt ~ c.n_soloid, data = dadosmisto1)
s   <- summary(fit)

r2  <- s$r.squared

coefs <- coef(s)
slope_row <- setdiff(rownames(coefs), "(Intercept)")
p  <- coefs[slope_row, "Pr(>|t|)"]

subtxt <- sprintf("R\u00B2 = %.2f; p %s",
                  r2,
                  ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

g <- ggplot(dadosmisto1, aes(c.n_soloid, log_produt)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, colour = "orange") +
  labs(
    x = "Soil C:N ratio",
    y = expression("log Biomass accumulation (g m"^-2*" year"^-1*")"),
    subtitle = subtxt
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.margin = unit(c(6,10,6,10), "mm"))

g

ggsave("~/01 Masters_LA/06 Figures/02 plots/logprodut_c.n_soloid.jpeg", g,
       width = 15, height = 10, units = "cm", dpi = 600, bg = "white")


## PCPS1

fit <- lm(log_produt ~ pcps1, data = dadosmisto1)
s   <- summary(fit)

r2  <- s$r.squared

coefs <- coef(s)
slope_row <- setdiff(rownames(coefs), "(Intercept)")
p  <- coefs[slope_row, "Pr(>|t|)"]

subtxt <- sprintf("R\u00B2 = %.2f; p= %s",
                  r2,
                  ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

g <- ggplot(dadosmisto1, aes(pcps1, log_produt)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, colour = "blue") +
  labs(
    x = "PCPS Axis 1",
    y = expression("log Biomass accumulation (g m"^-2*" year"^-1*")"),
    subtitle = subtxt
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.margin = unit(c(6,10,6,10), "mm"))

g

ggsave("~/01 Masters_LA/06 Figures/02 plots/logprodut_pcps1.jpeg", g,
       width = 15, height = 10, units = "cm", dpi = 600, bg = "white")


### FACET WRAP ###

#SR

g <- dadosmisto1 %>%
  ggplot(aes(sr, log_produt)) +
  geom_point(size = 3, alpha = 0.5) + 
  geom_smooth(method = "lm", se = FALSE) + 
  facet_wrap(~site) +
  labs(x = "Species richness (SR)",
       y = expression("Biomass accumulation (g m"^-2*" year"^-1*")")) +
  theme_bw(base_size = 11)

ggsave("~/01 Masters_LA/06 Figures/02 plots/facetwrap_SR.jpeg",
       plot = g,
       width = 13, height = 10, units = "cm", dpi = 600, bg = "white")


############################
### Plot Partial Effects ###
############################

# -------
# Model
# -------
m12 <- lmer(log_produt ~ c.n_soloid + SR + pcps1ab + (1 | site),
            data = dadosmisto1, REML = FALSE)

# p-values
pv <- summary(m12)$coefficients[, "Pr(>|t|)"]
p_txt_en <- sprintf("p(C:N)=%.3f | p(sr)=%.3f | p(PCPS1)=%.3f",
                    pv["c.n_soloid"], pv["sr"], pv["pcps1"])

# ------------------------------
# Data frame with english names
# ------------------------------
dados_plot <- dadosmisto1 %>%
  rename(
    soil_C_N_ratio      = c.n_soloid,
    species_richness_SR = sr,
    PCPS1               = pcps1,
    log_productivity    = log_produt
  )

# --------------------------------------------------------------------
# Marginal Effects (only fixed)
# --------------------------------------------------------------------
eff_c   <- ggpredict(m12, terms = "c.n_soloid", type = "fixed")
eff_sr  <- ggpredict(m12, terms = "sr",         type = "fixed")
eff_pc1 <- ggpredict(m12, terms = "pcps1",    type = "fixed")

# --------------------------------------------------------------------
# Plot with points + line + CI
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

# Individual plots
g_cn   <- plot_orig_pts(eff_c,   dados_plot, "soil_C_N_ratio",      "Soil C:N ratio")
g_sr   <- plot_orig_pts(eff_sr,  dados_plot, "species_richness_SR", "Species Richness (SR)") +
  ylab(NULL)  # remove y-label from the right painel
g_pcps <- plot_orig_pts(eff_pc1, dados_plot, "PCPS1",               "First PCPS axis")

# design: A (top-left), B (top-right); C (bottom-left), empty right
design <- "
AB
C.
"

#
fig_or <- wrap_plots(
  A = g_cn,        
  B = g_sr,        
  C = g_pcps,      
  design = design,
  heights = c(1, 0.55)   
) + plot_annotation(
  title = "Partial effects with observed points — original scale",
  subtitle = p_txt_en,
  tag_levels = "a",        
  tag_prefix = "(", tag_suffix = ")"
)

fig_or
ggsave("~/01 Masters_LA/06 Figures/02 plots/fig_m12_partial_effects_original_EN.png", fig_or, width = 11, height = 6.4, dpi = 300)

######################################
######### PCPS - without CSA #########
######################################

# input the sample species list without CSA
example <- read.csv("01 Datasets/01_raw_data/filogenia_total_Sem_CSA.csv")

### generate a phylogeny for the sample species list
tree <- phylo.maker(example, tree = GBOTB.extended.TPL,output.sp.list = TRUE,nodes = nodes.info.1.TPL, scenarios="S3")
tree_ok <- tree$scenario.3
summary(tree_ok)

phylo.dissim = cophenetic(tree_ok)
std.pdis = (phylo.dissim-min(phylo.dissim))/(max(phylo.dissim)-min(phylo.dissim))
pdis.ord = std.pdis[order(row.names(std.pdis)), order(row.names(std.pdis))]

# install.packages("PCPS")

### Importing the family groups
# Define Phylogenetic groups
grupos <- example$family

# Importing the community matrix
comunidade <- read.csv("01 Datasets/01_raw_data/comunidade_abund_sem_csa.csv", row.names = 1,header = T, sep = ";")

# Running PCPS
res <- pcps(comunidade, pdis.ord)
summary(res)

# Calculating explained variance
eig <- res$values$Eigenvalue
var_exp <- eig / sum(eig)
pcps1_var <- round(var_exp[1] * 100, 2)
pcps2_var <- round(var_exp[2] * 100, 2)

# Scores for sites and species
scores_pcps <- scores.pcps(res, choice = c(1,2))
scores_sites <- scores_pcps$scores.sites
scores_species <- scores_pcps$scores.species

# Group arrows (families)
apg_arrows <- apply(scores_species, 2, function(x) tapply(x, list(grupos), mean)) %>%
  as.data.frame()
apg_arrows$angle <- atan2(apg_arrows$pcps.1, apg_arrows$pcps.2) * (180/pi)
apg_arrows$labels <- c("Anac", "Anno", "Apocy", "Aster", "Bigno", "Borag", "Canna", 
                       "Caric", "Euphorb", "Faba", "Lami", "Laura", "Lecyt", "Malva", 
                       "Melast", "Melia", "Myrta", "Petive","Primula", "Rubi", 
                       "Ruta", "Sapind", "Solan", "Urtic", "Verbe")

# Merging the data into a single data.frame
df_plot <- cbind(as.data.frame(scores_sites),
                 produt = dadosmisto1$produt_z_g.ano,
                 SESPD = dadosmisto1$sespd,
                 SR = dadosmisto1$sr)

## only sespd
ggplot(data = df_plot, aes(x = pcps.1, y = pcps.2)) +
  geom_point(aes(color = SESPD), size = 3) +
  scale_colour_gradientn(
    colors = c("red", "pink", "#FAFD77", "black"),
    values = scales::rescale(c(-4, -2, 0, 2)),
    limits = c(-5, 2),
    oob = scales::squish,
    name = "sesPD"
  ) +
  labs(
    x = paste0("PCPS1 (", pcps1_var, "%)"),
    y = paste0("PCPS2 (", pcps2_var, "%)"),
    title = "Phylogenetic composition of the study areas"
  ) +
  geom_label(
    data = apg_arrows,
    aes(x = pcps.1, y = pcps.2, label = labels),
    fontface = "bold", fill = "grey", alpha = 0.5,
    color = "black", size = 5
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 13),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

ggsave("~/01 Masters_LA/06 Figures/02 plots/pcps_SESPD_semCSA.jpeg", width = 15, height = 10, dpi = 300, units = "in")

### only biomass

ggplot(data = df_plot, aes(x = pcps.1, y = pcps.2)) +
  geom_point(aes(color = produt), size = 4) +
  scale_colour_gradientn(
    colors = c("black", "lightyellow", "red"),
    name = "Biomass \nper year (g/y)"
  ) +
  labs(
    x = paste0("PCPS1 (", pcps1_var, "%)"),
    y = paste0("PCPS2 (", pcps2_var, "%)"),
    title = "Phylogenetic composition of the study areas"
  ) +
  geom_label(
    data = apg_arrows,
    aes(x = pcps.1, y = pcps.2, label = labels),
    fontface = "bold", fill = "grey", alpha = 0.5,
    color = "black", size = 5
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 13),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

ggsave("~/01 Masters_LA/06 Figures/02 plots/pcps_biomass_semCSA.jpeg", width = 15, height = 10, dpi = 300, units = "in")

################################################
## Family contribution (%) to PCPS1 and PCPS2 ##
################################################

example$species_us <- gsub(" ", "_", example$species)
species_in_scores <- rownames(scores_species)

family_vec <- example$family[
  match(species_in_scores, example$species_us)
]

fam_centroids <- aggregate(
  scores_species[, c("pcps.1", "pcps.2")],
  by = list(Family = family_vec),
  FUN = mean
)

fam_ss_pc1 <- tapply(scores_species[, "pcps.1"]^2,
                     family_vec,
                     sum)

fam_ss_pc2 <- tapply(scores_species[, "pcps.2"]^2,
                     family_vec,
                     sum)
fam_contrib <- tibble(
  Family = names(fam_ss_pc1),
  SS_PC1 = as.numeric(fam_ss_pc1),
  Pct_PC1 = 100 * SS_PC1 / sum(SS_PC1),
  SS_PC2 = as.numeric(fam_ss_pc2[names(fam_ss_pc1)]),
  Pct_PC2 = 100 * SS_PC2 / sum(SS_PC2)
)

fam_summary <- fam_centroids %>%
  left_join(fam_contrib, by = "Family") %>%
  arrange(desc(Pct_PC1)) 
head(fam_summary, 10)
# Fabaceae #1

###########################
### PHYLOGENETIC SIGNAL ###
###########################

# Without CSA

dados <- read.csv("01 Datasets/01_raw_data/prod_spp_semcsa.csv", row.names = 1,header = T, sep = ";")

tree_teste <- multi2di(tree_ok)

resultado <- multiPhylosignal(dados, tree_teste) # k = 0.1805496 FRACO SINAL

## DAP ### WITH CSA ## ALTER THIS LATER

dados <- read.csv("01 Datasets/01_raw_data/dap_spp.csv", row.names = 1,header = T, sep = ";")

tree_teste <- multi2di(tree_ok)

resultado <- multiPhylosignal(dados, tree_teste) # k = 0.2499149 FRACO SINAL

#######################################
########## PCA GRANULOMETRY ###########
#######################################

granulometria <- read_excel("01 Datasets/01_raw_data/granulometria.xlsx")

dados_granulometria <- granulometria[, -1]

pca_resultado <- prcomp(dados_granulometria, scale. = TRUE)
summary(pca_resultado)

screeplot(pca_resultado, main = "Screeplot - PCA")

escores_pca <- pca_resultado$x
head(escores_pca)

granulometria_pca <- cbind(granulometria, escores_pca)

granulometria_pca <- granulometria_pca %>%
  select(site, PC1, PC2)


# Plot (PC1 vs PC2)
ggplot(granulometria_pca, aes(x = PC1, y = PC2, color = site)) +
  geom_point(size = 4) +
  labs(title = "PCA: PC1 vs PC2", x = "Componente Principal 1 (PC1)", y = "Componente Principal 2 (PC2)") +
  theme_minimal() +
  scale_color_discrete(name = "Site") +
  theme(plot.title = element_text(hjust = 0.5))

# Variáveis de granulometria (para as setas)
variaveis_granulometria <- data.frame(
  x = pca_resultado$rotation[, 1],
  y = pca_resultado$rotation[, 2],
  nome = colnames(dados_granulometria)
)

# Plotar com ggplot2
ggplot(granulometria_pca, aes(x = PC1, y = PC2)) +
  # Adicionar os pontos com nome do site
  geom_point(size = 4) +
  geom_text(aes(label = site), vjust = -1.5, hjust = 0.5) +
  
  # Adicionar setas para as variáveis de granulometria
  geom_segment(data = variaveis_granulometria, aes(x = 0, y = 0, xend = x, yend = y), 
               arrow = arrow(type = "closed", length = unit(0.2, "inches")), 
               color = "red", size = 1) +
  geom_text(data = variaveis_granulometria, aes(x = x, y = y, label = nome), 
            vjust = -1, color = "red") +
  
  # Labels e tema
  labs(title = "PCA: PC1 vs PC2", x = "Componente Principal 1 (PC1)", y = "Componente Principal 2 (PC2)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("~/01 Masters_LA/00 MASTERS-DATA/01 Datasets/01_raw_data/PCA_granulometria.jpeg", width = 10, height = 8, dpi = 300)

cor(dados_granulometria, pca_resultado$x)
# optional
write_xlsx(granulometria_pca, "~/01 Masters_LA/00 MASTERS-DATA/01 Datasets/02_processed_data/granulometria_pca.xlsx")


########################################
########### PCA - NUTRIENTS ############
########################################

pca_nutrientes <- read_excel("~/01 Masters_LA/00 MASTERS-DATA/01 Datasets/01_raw_data/pca_nutrientes.xlsx")

dados_nutrientes <- pca_nutrientes[, -1]

pca_resultado <- prcomp(dados_nutrientes, scale. = TRUE)
summary(pca_resultado)

screeplot(pca_resultado, main = "Screeplot - PCA")

escores_pca <- pca_resultado$x
head(escores_pca)

nutrientes_pca <- cbind(pca_nutrientes, escores_pca)

nutrientes_pca <- nutrientes_pca %>%
  select(site, PC1, PC2)

# Create a data frame with the variable loadings in the first two components
variaveis_nutrientes <- as.data.frame(pca_resultado$rotation[, 1:2])
variaveis_nutrientes$nome <- rownames(variaveis_nutrientes) 
colnames(variaveis_nutrientes) <- c("x", "y", "nome") # Renomear colunas

variaveis_nutrientes <- variaveis_nutrientes %>%
  mutate(
    hjust_pos = ifelse(x > 0, 0.2, -0.5),
    vjust_pos = ifelse(y > 0, -0.8, 1)
  )

# PLOT
ggplot(nutrientes_pca, aes(x = PC1, y = PC2)) +
  geom_point(size = 4) +
  geom_text(aes(label = site), vjust = -1.5, hjust = 0.5) +
  geom_segment(data = variaveis_nutrientes, aes(x = 0, y = 0, xend = x, yend = y),
               arrow = arrow(type = "closed", length = unit(0.2, "inches")),
               color = "red", size = 1) +
  geom_text(data = variaveis_nutrientes, 
            aes(x = x, y = y, label = nome), 
            color = "red", size = 5,
            vjust = variaveis_nutrientes$vjust_pos, 
            hjust = variaveis_nutrientes$hjust_pos) +
  labs(title = "PCA: Nutrientes no solo",
       x = "Componente Principal 1 (PC1)",
       y = "Componente Principal 2 (PC2)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("~/01 Masters_LA/06 Figures/02 plots/PCA_nutrientes.jpeg", width = 10, height = 8, dpi = 300)

# Calculate the correlation between the original variables and the PCs.
correlacoes <- cor(dados_nutrientes, pca_resultado$x)
print(correlacoes)

# Salvar os resultados em Excel
write_xlsx(nutrientes_pca, "01 Datasets/02_processed_data/nutrientes_pca.xlsx")
write_xlsx(as.data.frame(correlacoes), "01 Datasets/02_processed_data/correlacoes_pca.xlsx")

##################################################
# SEM - PCA Climatic Variables (ppt, tmax, tmin) #
##################################################

# 2. Carregar os dados (exemplo com dados fictícios)

dados <- read.csv("01 Datasets/01_raw_data/dadosmisto.csv",
                       header = TRUE,
                       sep = ";")

dados1 <- dados[-c(1:7), ]

# scale
dados_padronizados <- dados %>% 
  select(ppt, tmax, tmin) %>% 
  scale() %>% 
  as.data.frame()

# PCA
pca_resultado <- PCA(dados_padronizados, graph = FALSE)

summary(pca_resultado) 
fviz_eig(pca_resultado) # Scree plot

# Extract cargas (loadings)


jpeg("~/01 Masters_LA/06 Figures/02 plots/pca_aridez.jpeg", width = 2000, height = 1600, res = 300)
print(fviz_pca_var(pca_resultado, 
                   col.var = "contrib",
                   gradient.cols = c("blue", "yellow", "red"),
                   repel = TRUE,
                   title = "Dryness Gradient"))
dev.off()

# PC1 extraction to use in SEM
dados$PC1_clima <- pca_resultado$ind$coord[, 1]

# Check the direction of PC1 (important for interpretation)
cor(dados_padronizados, dados$PC1_clima)

# 9. Salvar os dados com o PC1
write.csv(dados, "01 Datasets/02_processed_data/dados_pc1aridez.csv", row.names = FALSE)

########################################
# Phylogenetic Analysis - HILL NUMBERS #
########################################

# install.packages(c("picante", "PhyloMeasures", "hillR", "ape", "phytools", "dplyr"))
#install.packages("hillR")

comunidade <- read.csv("01 Datasets/01_raw_data/comunidade_abund_sem_csa.csv", row.names = 1,header = T, sep = ";")

# a) hillR
hill_metrics <- data.frame(
  parcela = rownames(comunidade),
  PD_q0 = hill_phylo(comunidade, tree_ok, q = 0),
  PD_q1 = hill_phylo(comunidade, tree_ok, q = 1),
  PD_q2 = hill_phylo(comunidade, tree_ok, q = 2)
)

write_xlsx(hill_metrics, "~/01 Masters_LA/00 MASTERS-DATA/01 Datasets/02_processed_data/hill_metrics.xlsx")


# LMMs

modelo_nulo <- lmer(log_produt ~ (1 | site), data = dadosmisto, REML = FALSE)
summary(modelo_nulo)
AICc(modelo_nulo)


misto1 <- lmer(log_produt ~ PD_q0 + (1 | site), data = dadosmisto1, REML = FALSE)
summary(misto1) # p-value =  0.011 *
r.squaredGLMM(misto1) # 0.6126289
AICc(misto1) # 75.58639

misto2 <- lmer(log_produt ~ PD_q1 + (1 | site), data = dadosmisto1, REML = FALSE)
summary(misto2) # p-value = 0.33
r.squaredGLMM(misto2) # 0.5093365
AICc(misto2) # 80.89568

misto3 <- lmer(log_produt ~ PD_q2 + (1 | site), data = dadosmisto1, REML = FALSE)
summary(misto3) # p-value = 0.24
r.squaredGLMM(misto3) # 0.5127555
AICc(misto3) # 80.4419

misto4 <- lmer(log_produt ~ PD_q0 + sespd + (1 | site), data = dadosmisto1, REML = FALSE)
summary(misto4) 
r.squaredGLMM(misto4)
AICc(misto4) 

