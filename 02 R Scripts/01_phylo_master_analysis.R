########################################
### SCRIPT FOR PHYLOGENETIC ANALYSIS ###
########################################

# last update: 2025/10/27 (YYYY/MM/DD)
# author: Laíla Arnauth

################################
#### PHYLOGENETIC DIVERSITY ####
################################

# Load all packages

library(readr)
library("V.PhyloMaker2")
library(writexl)
library(picante)
library(ape)
library(vegan)
library(readxl)
library(tidyverse)
library(corrplot)
library(lmerTest) #MODELO MISTO
library(MuMIn) #AICc
library(car)

# Install V.PhyloMaker2
# install.packages("devtools")
# devtools::install_github("jinyizju/V.PhyloMaker2")

# input the sample species list
example <- read_csv("02 Datasets/01 raw_data/filogenia_total.csv")
View(example)

### generate a phylogeny for the sample species list
tree <- phylo.maker(example, tree = GBOTB.extended.TPL,output.sp.list = TRUE,nodes = nodes.info.1.TPL, scenarios="S3")
tree_ok <- tree$scenario.3
summary(tree_ok)

# Plot Phylo Tree

plot.phylo(tree_ok, type = "tidy", show.tip.label = TRUE, show.node.label = FALSE,edge.color = "gray",edge.width = 1,tip.color = "black")

png("arvore_filog.png",width = 12, height = 10, units = "in", res = 300)
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

comunidade <- read.csv("02 Datasets/01 raw_data/comunidade_abund.csv",
  row.names = 1,
  header = TRUE,
  sep = ";")

### Finally, phylogenetic diversity is calculated w/ picante
colnames(comunidade)==colnames(pdis.ord) # checking if they match

# PD FAITH

pd.faith <- pd(comunidade,tree_ok,include.root = FALSE)
view(pd.faith)
write_xlsx(pd.faith, "~/01 Masters_LA/MASTERS-THESIS/02 Datasets/02 processed_data/PD_Faith.xlsx")

pd_only <- pd.faith %>%
  tibble::rownames_to_column(var = "site") %>%
  dplyr::select(-SR)


# MPD

mpd_df <- tibble::tibble(
  site = rownames(comunidade),
  mpd  = as.numeric(mpd(comunidade, pdis.ord, abundance.weighted = TRUE))
)
write_xlsx(mpd_df, "~/01 Masters_LA/MASTERS-THESIS/02 Datasets/02 processed_data/MPD.xlsx")

# MNTD (species are more broadly or finely spaced in a community compared to another?)

mntd_df <- tibble::tibble(
  site = rownames(comunidade),
  mntd = as.numeric(mntd(comunidade, pdis.ord, abundance.weighted = TRUE))
)

write_xlsx(mntd_df, "~/01 Masters_LA/MASTERS-THESIS/02 Datasets/02 processed_data/MNTD.xlsx")

### Standardized effect size of PD ##### Runs = 99 # ape
# tree_clean <- drop.tip(tree_ok, tip = which(is.na(tree_ok$tip.label))) #OPTIONAL

ses_pd <- ses.pd(comunidade, tree_ok, null.model="richness", runs = 999,iterations = 1000,include.root = FALSE)
ses_pd <- as.data.frame(ses_pd)

write_xlsx(ses_pd, "~/01 Masters_LA/MASTERS-THESIS/02 Datasets/02 processed_data/ses_pd.xlsx")

ses_pd <- ses_pd |> 
  tibble::rownames_to_column(var = "site")

sespd_only <- ses_pd |> 
  transmute(site, sespd = pd.obs.z)
view(sespd_only)

### Standardized effect size of MPD ##### Runs = 999

ses_mpd <- ses.mpd(comunidade, pdis.ord,null.model="richness")
ses_mpd <- as.data.frame(ses_mpd)

write_xlsx(ses_mpd, "~/01 Masters_LA/MASTERS-THESIS/02 Datasets/02 processed_data/ses_mpd.xlsx")

ses_mpd <- ses_mpd |> 
  tibble::rownames_to_column(var = "site")

sesmpd_only <- ses_mpd |> 
  transmute(site, sesmpd = mpd.obs.z)
view(sesmpd_only)

### Standardized effect size of MNTD ##### Runs = 999

ses_mntd <- ses.mntd(comunidade, pdis.ord,null.model="richness")
ses_mntd <- as.data.frame(ses_mntd)
write_xlsx(ses_mntd, "~/01 Masters_LA/MASTERS-THESIS/02 Datasets/02 processed_data/ses_mntd.xlsx")

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
write_xlsx(PSR, "~/01 Masters_LA/MASTERS-THESIS/02 Datasets/02 processed_data/PSR.xlsx")

# Phylogenetic Species Evenness: is the metric PSV modified to incorporate relative species abundances
PSE <- pse(comunidade,tree_ok)
write_xlsx(PSE, "~/01 Masters_LA/MASTERS-THESIS/02 Datasets/02 processed_data/PSE.xlsx")

# Phylogenetic species clustering: is a metric of the branch tip clustering of species across a sample’s phylogeny

PSC <- psc(comunidade,tree_ok)
write_xlsx(PSC, "~/01 Masters_LA/MASTERS-THESIS/02 Datasets/02 processed_data/PSC.xlsx")

###########################
### TAXONOMIC DIVERSITY ###
###########################

comunidade_diversity <- read.csv("02 Datasets/01 raw_data/comunidade_diversity.csv",
                                 row.names = 1,
                                 header = TRUE,
                                 sep = ";")
# Simpson
unbias_simp_df <- tibble::tibble(
  site = rownames(comunidade_diversity),
  unbias_simp = as.numeric(simpson.unb(comunidade_diversity))
)
write_xlsx(unbias.simp.df, "~/01 Masters_LA/MASTERS-THESIS/02 Datasets/02 processed_data/Unbias_Simp.xlsx")

#Shannon
shannon_df <- tibble::tibble(
  site = rownames(comunidade_diversity),
  shannon = as.numeric(diversity(comunidade_diversity, index = "shannon"))
)

write_xlsx(shannon_df, "~/01 Masters_LA/MASTERS-THESIS/02 Datasets/02 processed_data/shannon.xlsx")

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
  "~/01 Masters_LA/MASTERS-THESIS/02 Datasets/02 processed_data/data_diversity.csv"
)

write_xlsx(data_diversity,
           "~/01 Masters_LA/MASTERS-THESIS/02 Datasets/02 processed_data/data_diversity.xlsx")

########################
### DATA EXPLORATION ###
########################

dadosmisto <- read.csv("02 Datasets/01 raw_data/dadosmisto.csv",
                                 header = TRUE,
                                 sep = ";")

dadosmisto1 <- dadosmisto[-c(1:7), ]

png("~/01 Masters_LA/MASTERS-THESIS/02 Datasets/04 exploratory_plots/hist_biomass.png", width = 800, height = 600)
hist(dadosmisto$biomassa_z_kg, main="Histograma Biomassa", xlab="Valores", col="blue", border="black")
dev.off()

png("~/01 Masters_LA/MASTERS-THESIS/02 Datasets/04 exploratory_plots/hist_logprodut.png", width = 800, height = 600)
hist(dadosmisto$log_produt, main="Histograma Log Produtividade", xlab="Valores", col="green", border="black")
dev.off()

# SHAPIRO TEST

shap_biomass <- shapiro.test(dadosmisto$biomassa_z_kg) # NOT NORMAL
dados_transformados <- log(dadosmisto$biomassa_z_kg + 1)  # Add 1 to avoid log(0)
dados_transformados.df <- as.data.frame(dados_transformados)
write_xlsx(dados_transformados.df, "~/01 Masters_LA/MASTERS-THESIS/02 Datasets/02 processed_data/logbiomassa.xlsx")
shapiro.test(dados_transformados) 

shap_produt <- shapiro.test(dadosmisto$produt_z_g.ano) # NÃO É NORMAL
dados_transformados <- log(dadosmisto$produt_z_g.ano + 1)  # Adicione 1 para evitar log(0)
dados_transformados.df <- as.data.frame(dados_transformados)
shapiro.test(dados_transformados) 
write_xlsx(dados_transformados.df, "~/01 Masters_LA/MASTERS-THESIS/02 Datasets/02 processed_data/logprodut.xlsx")


# DENSITY

plot(density(dadosmisto$log_biomass), main="Densidade dos Dados")
plot(density(dadosmisto$log_produt), main="Densidade dos Dados")

# CHECKING FOR OUTLIERS

png("~/01 Masters_LA/MASTERS-THESIS/02 Datasets/04 exploratory_plots/boxplot_LogBiomassa.png",
    width = 800, height = 600)
boxplot(dadosmisto$log_biomass, 
        main = "Boxplot - Log Biomass (kg)",
        col = "#8B008B")
dev.off()

png("~/01 Masters_LA/MASTERS-THESIS/02 Datasets/04 exploratory_plots/boxplot_LogProdutividade.png",
    width = 800, height = 600)
boxplot(dadosmisto$log_produt, 
        main = "Boxplot - Productivity (g/cm²/year)", 
        col = "#32CD32")
dev.off()

############# CORRELATION MATRIX - INDEPENDENT VARIABLES (X) #############

# tidyverse

dadosmisto <- read.csv("02 Datasets/01 raw_data/dadosmisto.csv",
                       header = TRUE,
                       sep = ";")
dadosmisto1 <- dadosmisto[-c(1:7), ] # WITHOUT CSA

# Defining the independent variables (X)

variaveis_x <- c("SR","unbias.simp","shannon","SESPDab", "SESMPDab", 
                 "SESMNTDab", "PSV", "PSR", "PSEab", "PSC","ppt", 
                 "tmax", "tmin", "pet", "vpd","mcwd","PC1_clima","altitude","declividade",
                 "silte","ph","valor_s","valor_t","PC1nutri","PC2nutri","pcps1ab","pcps2ab","faba","PD_q0")

variaveis_x <- c("SR","silte","pcps1ab","PD_q0","declividade")

variaveis_x <- c("SR","SESPDab", "SESMPDab", 
                 "SESMNTDab", "pcps1ab","pcps2ab","faba","PD_q0")

variaveis_x <- c("ppt","tmax", "tmin", "pet", "vpd","mcwd","PC1_clima","altitude","declividade",
                 "silte","ph","valor_s","valor_t","PC1nutri","PC2nutri")

# Creating a new dataframe containing only the X variables.
dados_x <- dadosmisto1 %>%
  select(all_of(variaveis_x))

# Calculating the correlation matrix with the dataframe dados_x
cor_matrix_x <- cor(dados_x, use = "complete.obs")
print(cor_matrix_x)

# Plot

png("~/01 Masters_LA/MASTERS-THESIS/02 Datasets/04 exploratory_plots/matriz.cor_environmental_withoutCSA.png", width = 1200, height = 800, res = 150)
corrplot(cor_matrix_x, method = "circle", type = "upper", 
         tl.col = "black", tl.cex = 0.8)
dev.off()



##########################
### LINEAR MIXED MODEL ###
##########################

dadosmisto <- read.csv("02 Datasets/01 raw_data/dadosmisto.csv",
                       header = TRUE,
                       sep = ";")

dadosmisto1 <- dadosmisto[-c(1:7), ]

# OFICIAL ### site as random variable

modelo_nulo <- lmer(log_produt ~ (1|site), data = dadosmisto1, REML = FALSE)

### SEM CSA ###

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

################
### GRÁFICOS ###
################

## ==========================================================
##  Coefplot do modelo m11 — estilo Restoration Ecology
## ==========================================================

library(ggplot2)
library(broom.mixed)
library(grid)  # para unit()


## 1) Padronizar as variáveis (média = 0, desvio padrão = 1)

dados_std <- dadosmisto1
dados_std$c.n_soloid <- scale(dadosmisto1$c.n_soloid)
dados_std$SR         <- scale(dadosmisto1$SR)
dados_std$silte      <- scale(dadosmisto1$silte)
dados_std$pcps1ab    <- scale(dadosmisto1$pcps1ab)


## 2) Rodar o modelo padronizado (igual ao m11 original)

m11_std <- lmer(log_produt ~ c.n_soloid + SR + pcps1ab + (1 | site),
                data = dados_std, REML = FALSE)

summary(m11_std)

## 3) Extrair coeficientes + IC95% do modelo padronizado
coef_df <- tidy(m11_std, effects = "fixed", conf.int = TRUE)

## 4) Remover intercepto e criar coluna de significância
coef_df_sem_intercepto <- subset(coef_df, term != "(Intercept)")
coef_df_sem_intercepto$significativo <- ifelse(coef_df_sem_intercepto$p.value < 0.05,
                                               "Significant", "Non-significant")

## 5) Rótulos e ordem desejada (SR topo → PCPS1 → C:N → Silt)
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
setwd("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/graficos")
ggsave("graf_coefplot_msel.jpg",
       plot = p, width = 14, height = 9, units = "cm",
       dpi = 600, bg = "white", limitsize = FALSE)



################
### PLOT GLM ###
################

library(writexl)
library(readxl)
library(tidyverse)
library(ggplot2)
setwd("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/graficos")
dadosmisto <- read_excel("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/dadosmisto.xlsx")


## SR 

fit <- lm(log_produt ~ SR, data = dadosmisto1)
summary(fit)

r2  <- summary(fit)$r.squared
p   <- coef(summary(fit))["SR", "Pr(>|t|)"]

subtxt <- sprintf("R\u00B2 = %.2f; p = %s",
                  r2,
                  ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))


g <- ggplot(dadosmisto1, aes(SR, log_produt)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, colour = "yellow") +
  labs(
    x = "Species richness (SR)",
    y = expression("log Biomass accumulation (g m"^-2*" year"^-1*")"),
    subtitle = subtxt
  ) +
  theme_minimal(base_size = 11)
g

ggsave("logprodut_SR.jpeg", g, width = 15, height = 10, units = "cm", dpi = 600)



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

ggsave("logprodut_silt.jpeg", g,
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

ggsave("logprodut_c.n_soloid.jpeg", g,
       width = 15, height = 10, units = "cm", dpi = 600, bg = "white")


## PCPS1

fit <- lm(log_produt ~ pcps1ab, data = dadosmisto1)
s   <- summary(fit)

r2  <- s$r.squared

coefs <- coef(s)
slope_row <- setdiff(rownames(coefs), "(Intercept)")
p  <- coefs[slope_row, "Pr(>|t|)"]

subtxt <- sprintf("R\u00B2 = %.2f; p= %s",
                  r2,
                  ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

g <- ggplot(dadosmisto1, aes(pcps1ab, log_produt)) +
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

ggsave("logprodut_pcps1.jpeg", g,
       width = 15, height = 10, units = "cm", dpi = 600, bg = "white")


### FACET WRAP ###
library(tidyverse)


#SR

g <- dadosmisto1 %>%
  ggplot(aes(SR, log_produt)) +
  geom_point(size = 3, alpha = 0.5) + 
  geom_smooth(method = "lm", se = FALSE) + 
  facet_wrap(~site) +
  labs(x = "Species richness (SR)",
       y = expression("Biomass accumulation (g m"^-2*" year"^-1*")")) +
  theme_bw(base_size = 11)

ggsave("facetwrap_SR.jpeg",
       plot = g,
       width = 13, height = 10, units = "cm", dpi = 600, bg = "white")

# PCPS1

g <- dadosmisto1 %>%
  ggplot(aes(pcps1ab, log_produt)) +
  geom_point(size = 3, alpha = 0.5) + 
  geom_smooth(method = "lm", se = FALSE) + 
  facet_wrap(~site) +
  labs(x = "PCPS Axis 1",
       y = expression("Biomass accumulation (g m"^-2*" year"^-1*")")) +
  theme_bw(base_size = 11)

ggsave("facetwrap_pcps1.jpeg",
       plot = g,
       width = 13, height = 10, units = "cm", dpi = 600, bg = "white")



##################################
### GRÁFICOS MODELOS MÚLTIPLOS ###
##################################

install.packages("sjPlot")
library(ggplot2)
library(sjPlot)

# Gerando gráficos de efeitos parciais para o modelo
plot_model(M_PSE_ppt_decliv_shan, type = "eff", terms = c("PSE", "ppt", "declividade", "shannon"))

# Gerando gráficos de efeitos parciais para o modelo
plot_model(M_PSE_ppt_decliv, type = "eff", terms = c("PSE", "ppt", "declividade"))


#############################
### GRÁFICOS MODELO MISTO ###
#############################

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

### PREDICT ####

# Carregue os pacotes
library(effects)
library(ggeffects)

# Gráfico de efeitos preditivos para PSE e SR
eff <- effect("silt*SR", misto1)
plot(eff)

library(ggplot2)

# Gerar o gráfico com ggpredict
pred <- ggpredict(misto1, terms = c("SR", "silt"))
plot(pred)

# Criar o gráfico
grafico <- plot(pred) +
  labs(title= "Predicted values of Productivity",x = "SR", y = "log Productivity (g/m²/ano)") +
  theme_minimal()

plot(grafico)

###
# Salvar o gráfico como PNG
setwd("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/graficos")
ggsave("Predict_SR_silte.png", plot = grafico, width = 15, height = 10, units = "cm",dpi = 300,bg = "white")

#### SR e PCPS

library(ggplot2)
library(ggeffects)

# Gerar previsões para a interação
pred_interacao <- ggpredict(model, terms = c("SR", "pcps1ab [-1, 0, 1]"))

# Plotar
ggplot(pred_interacao, aes(x = x, y = predicted, color = group)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2) +
  labs(x = "Riqueza de Espécies (SR)", 
       y = "Log(Produtividade)", 
       color = "PCPS1",
       title = "Efeito da riqueza sobre a produtividade em diferentes contextos filogenéticos") +
  theme_bw()

ggpredict(model, terms = c("pcps1ab", "SR [5, 10, 15]")) |> plot()

library(ggeffects)
library(ggplot2)

# Gerar previsões para SR (5 a 25) e pcps1ab (baixo, médio, alto)
pred_data <- ggpredict(model, terms = c("SR [5:25 by=5]", "pcps1ab [-0.5, 0, 0.5]"))

# Plotar
ggplot(pred_data, aes(x = x, y = predicted, color = group)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2) +
  labs(
    x = "Riqueza de Espécies (SR)", 
    y = "Log(Produtividade)", 
    color = "PCPS1",
    fill = "PCPS1",
    title = "Efeito não-linear de SR e pcps1ab na produtividade"
  ) +
  theme_bw()

######################################
######### PCPS - without CSA #########
######################################

# load the packages

library("V.PhyloMaker2")
library(writexl)

# input the sample species list
example <- read.csv("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/filogenia_total_Sem_CSA.csv")

### generate a phylogeny for the sample species list
tree <- phylo.maker(example, tree = GBOTB.extended.TPL,output.sp.list = TRUE,nodes = nodes.info.1.TPL, scenarios="S3")
tree_ok <- tree$scenario.3
summary(tree_ok)

# A função 'cophenetic' calcula a dissimilaridade com base nas distâncias filogenéticas da árvore
phylo.dissim = cophenetic(tree_ok)
# Agora padronizamos as distâncias em uma escala de variação de 0 a 1
std.pdis = (phylo.dissim-min(phylo.dissim))/(max(phylo.dissim)-min(phylo.dissim))
# Ordenamos as espécies em ordem alfabética, como na matriz da comunidade
pdis.ord = std.pdis[order(row.names(std.pdis)), order(row.names(std.pdis))]

install.packages("PCPS")
library(PCPS)
library(ggplot2)
library(tidyverse)
setwd("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/graficos")

### Importing the family groups
# Define Phylogenetic groups
grupos <- example$family

# Importing the community matrix
comunidade <- read.csv("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/comunidade_abund_sem_csa.csv", row.names = 1,header = T, sep = ";")

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
                 SESPDab = dadosmisto1$SESPDab,
                 SR = dadosmisto1$SR)

# Now the ggplot works

ggplot(data = df_plot, aes(x = pcps.1, y = pcps.2)) +
  geom_point(aes(color = SESPDab, size = produt)) +
  scale_colour_gradientn(
    colors = c("red", "lightyellow", "black"),  # azul → amarelo → vermelho
    name = "sesPD"
  ) +
  scale_size_continuous(
    name = "Biomass 
per year (g/y)", 
    range = c(1, 10),  # controla visualmente o menor e maior tamanho
    breaks = scales::pretty_breaks(n = 4)  # define 4 categorias de tamanho
  ) +
  labs(x = paste0("PCPS1 (", pcps1_var, "%)"), 
       y = paste0("PCPS2 (", pcps2_var, "%)"),
       title = "Phylogenetic composition of the study areas (w/ CSA)") +
  geom_text(data = apg_arrows,
            aes(x = pcps.1, y = pcps.2, label = labels), 
            fontface = "bold", color = "black", size = 1) +
  geom_label(data = apg_arrows,
             aes(x = pcps.1, y = pcps.2, label = labels), 
             fontface = "bold", fill = "grey", alpha = 0.5,
             color = "black", size = 5) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 13),
    legend.title = element_text(size = 14),   # título da legenda maior
    legend.text = element_text(size = 12)     # texto das categorias maior
  )

ggsave("pcps_semCSA.jpeg", width = 15, height = 10, dpi = 300, units = "in")

## only sespd
ggplot(data = df_plot, aes(x = pcps.1, y = pcps.2)) +
  geom_point(aes(color = SESPDab), size = 3) +
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

ggsave("pcps_biomass_semCSA.jpeg", width = 15, height = 10, dpi = 300, units = "in")

########################################
### SINAL FILOGENÉTICO PRODUTIVIDADE ###
########################################

## BIOMASSA ### SERÁ QUE BIOMASSA TEM UM SINAL FILOGENÉTICO?
# Média de biomassa por espécie

dados <- read.csv("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/prod_spp_semcsa.csv", row.names = 1,header = T, sep = ";")

library(ape)

tree_teste <- multi2di(tree_ok)

library(picante)

resultado <- multiPhylosignal(dados, tree_teste) # k = 0.1805496 FRACO SINAL

## DAP ### SERÁ QUE DAP TEM UM SINAL FILOGENÉTICO?
# Média de DAP por espécie

dados <- read.csv("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/dap_spp.csv", row.names = 1,header = T, sep = ";")

library(ape)

tree_teste <- multi2di(tree_ok)

library(picante)

resultado <- multiPhylosignal(dados, tree_teste) # k = 0.2499149 FRACO SINAL

###########################################
########## PCA DA GRANULOMETRIA ###########
###########################################

library(ggplot2)
library(readxl)
library(tidyverse)
setwd("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/graficos")

# Ler os dados
granulometria <- read_excel("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/granulometria.xlsx")

# Remover a coluna 'site' para fazer o PCA com apenas as variáveis numéricas
dados_granulometria <- granulometria[, -1]

# Realizar o PCA usando a função base prcomp()
pca_resultado <- prcomp(dados_granulometria, scale. = TRUE)
summary(pca_resultado)

# Plotar o Screeplot para ver a variância explicada por cada componente
screeplot(pca_resultado, main = "Screeplot - PCA")

# Obter os escores dos componentes principais (PCs)
escores_pca <- pca_resultado$x
head(escores_pca)

# Adicionar os escores PCA ao dataframe original
granulometria_pca <- cbind(granulometria, escores_pca)

# Selecionar os eixos mais explicativos (exemplo: PC1 e PC2)
granulometria_pca <- granulometria_pca %>%
  select(site, PC1, PC2)


# Plotar os dois primeiros componentes principais (PC1 vs PC2) com ggplot2
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

ggsave("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/graficos/PCA_granulometria.jpeg", width = 10, height = 8, dpi = 300)

cor(dados_granulometria, pca_resultado$x)


library(writexl)

write_xlsx(granulometria_pca, "C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/granulometria_pca.xlsx")


###########################################
########### PCA DOS NUTRIENTES ############
###########################################

library(ggplot2)
library(readxl)
library(tidyverse)
library(writexl)
setwd("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/graficos")

# Ler os dados
pca_nutrientes <- read_excel("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/pca_nutrientes.xlsx")

# Remover a coluna 'site' para fazer o PCA com apenas as variáveis numéricas
dados_nutrientes <- pca_nutrientes[, -1]

# Realizar o PCA usando a função base prcomp()
pca_resultado <- prcomp(dados_nutrientes, scale. = TRUE)
summary(pca_resultado)

# Plotar o Screeplot para ver a variância explicada por cada componente
screeplot(pca_resultado, main = "Screeplot - PCA")

# Obter os escores dos componentes principais (PCs)
escores_pca <- pca_resultado$x
head(escores_pca)

# Adicionar os escores PCA ao dataframe original
nutrientes_pca <- cbind(pca_nutrientes, escores_pca)

# Selecionar os eixos mais explicativos (exemplo: PC1 e PC2)
nutrientes_pca <- nutrientes_pca %>%
  select(site, PC1, PC2)

# Criar data frame com os loadings das variáveis nos dois primeiros componentes
variaveis_nutrientes <- as.data.frame(pca_resultado$rotation[, 1:2])
variaveis_nutrientes$nome <- rownames(variaveis_nutrientes) # Adicionar nome das variáveis
colnames(variaveis_nutrientes) <- c("x", "y", "nome") # Renomear colunas

# Ajustar a posição dos rótulos dependendo da direção do eixo
variaveis_nutrientes <- variaveis_nutrientes %>%
  mutate(
    hjust_pos = ifelse(x > 0, 0.2, -0.5),  # Deslocar para a direita se PC1 positivo, para a esquerda se negativo
    vjust_pos = ifelse(y > 0, -0.8, 1)  # Ajuste vertical baseado no valor de PC2
  )

# Gráfico PCA ajustado
ggplot(nutrientes_pca, aes(x = PC1, y = PC2)) +
  # Pontos representando os sites
  geom_point(size = 4) +
  geom_text(aes(label = site), vjust = -1.5, hjust = 0.5) +
  
  # Adicionar setas das variáveis
  geom_segment(data = variaveis_nutrientes, aes(x = 0, y = 0, xend = x, yend = y),
               arrow = arrow(type = "closed", length = unit(0.2, "inches")),
               color = "red", size = 1) +
  
  # Ajustar rótulos das variáveis com deslocamento
  geom_text(data = variaveis_nutrientes, 
            aes(x = x, y = y, label = nome), 
            color = "red", size = 5,
            vjust = variaveis_nutrientes$vjust_pos, 
            hjust = variaveis_nutrientes$hjust_pos) + # Deslocamento para a esquerda ou direita
  
  # Labels e tema
  labs(title = "PCA: Nutrientes no solo",
       x = "Componente Principal 1 (PC1)",
       y = "Componente Principal 2 (PC2)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Salvar o gráfico
ggsave("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/graficos/PCA_nutrientes.jpeg", width = 10, height = 8, dpi = 300)

# Calcular correlação entre variáveis originais e os PCs
correlacoes <- cor(dados_nutrientes, pca_resultado$x)
print(correlacoes)

# Salvar os resultados em Excel
write_xlsx(nutrientes_pca, "C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/nutrientes_pca.xlsx")
write_xlsx(as.data.frame(correlacoes), "C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/correlacoes_pca.xlsx")

# ------------------------------------------
# PCA para variáveis climáticas (ppt, tmax, tmin)
# Script para uso em SEM - Laíla Iglesias
# ------------------------------------------

# 1. Carregar pacotes
library(tidyverse)   # Manipulação de dados e visualização
library(FactoMineR)  # PCA avançada
library(factoextra)  # Visualização da PCA

# 2. Carregar os dados (exemplo com dados fictícios)

dados <- read_excel("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/dadosmisto.xlsx")

# 3. Padronizar as variáveis (média = 0, DP = 1)
dados_padronizados <- dados %>% 
  select(ppt, tmax, tmin) %>% 
  scale() %>% 
  as.data.frame()

# 4. Executar a PCA
pca_resultado <- PCA(dados_padronizados, graph = FALSE)

# 5. Visualizar resultados
summary(pca_resultado)  # Proporção de variância explicada
fviz_eig(pca_resultado) # Scree plot

# 6. Extrair cargas (loadings) das variáveis

setwd("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/graficos")
jpeg("pca_aridez.jpeg", width = 2000, height = 1600, res = 300)
print(fviz_pca_var(pca_resultado, 
                   col.var = "contrib",
                   gradient.cols = c("blue", "yellow", "red"),
                   repel = TRUE,
                   title = "Dryness Gradient"))
dev.off()

# 7. Obter o PC1 para uso no SEM
dados$PC1_clima <- pca_resultado$ind$coord[, 1]  # Extrai scores do PC1

# 8. Verificar a direção do PC1 (importante para interpretação)
cor(dados_padronizados, dados$PC1_clima)

# 9. Salvar os dados com o PC1
setwd("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/TRABALHO EXCEL")
write.csv(dados, "dados_pc1aridez.csv", row.names = FALSE)

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



###############
### COM CSA ###
###############
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



##########################################
############ PREDICT EFFECTS #############
##########################################

library(readxl)
library(car)
library(ggeffects)
library(ggplot2)
library(lme4)

dadosmisto <- read_excel("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/dadosmisto.xlsx")

misto.Produt.PSE.ppt.SR.Site <- lmer(log_produt ~ PSE + ppt + SR + (1 | site), data = dadosmisto, REML = FALSE)

# Plotar resíduos do modelo
plot(misto.Produt.PSE.ppt.SR.Site)

# Gerar efeitos marginais para a variável PSE
effects_PSE <- ggpredict(misto.Produt.PSE.ppt.SR.Site, terms = "PSE")
effects_SR <- ggpredict(misto.Produt.PSE.ppt.SR.Site, terms = "SR")
effects_ppt <- ggpredict(misto.Produt.PSE.ppt.SR.Site, terms = "ppt")

# Plotar os efeitos marginais com alterações no título e nos eixos
grafico <- plot(effects_SR) 
setwd("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/graficos")
ggsave("grafico_Predicted.PSE.SR.jpeg", plot = grafico, device = "jpeg", width = 8, height = 6, dpi = 300)

# PSE
grafico <- plot(effects_PSE) +
  ggtitle("Predicted values of productivity (best LMM model)") +
  xlab("PSE") +
  ylab("Productivity (g/m²/year)")

ggsave("Predicted_PSE.jpeg", plot = grafico, device = "jpeg", width = 8, height = 6, dpi = 300)

##### SR
grafico <- plot(effects_SR) +
  ggtitle("Predicted values of productivity (best LMM model)") +
  xlab("SR") +
  ylab("Productivity (g/m²/year)")


ggsave("Predicted_SR.jpeg", plot = grafico, device = "jpeg", width = 8, height = 6, dpi = 300)


##### SR
grafico <- plot(effects_ppt) +
  ggtitle("Predicted values of productivity (best LMM model)") +
  xlab("ppt") +
  ylab("Productivity (g/m²/year)")


ggsave("Predicted_ppt.jpeg", plot = grafico, device = "jpeg", width = 8, height = 6, dpi = 300)


# Criar gráficos de resíduos parciais para as variáveis do modelo

jpeg("crPlots_M6_ppt_dec.jpeg", width = 800, height = 600)

crPlots(M6.ppt.dec)

dev.off()

####################################
####### STRUCTURAL DIVERSITY #######
####################################

# Pacote necessário para manipulação de dados
if (!require("dplyr")) install.packages("dplyr")
library(dplyr)
library(readxl)

data <- read_excel("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/TRABALHO EXCEL/dap.xlsx")

# Função para calcular o coeficiente de Gini
gini_coefficient <- function(data) {
  data <- sort(data)  # Ordenar os dados
  n <- length(data)  # Número de observações
  mean_value <- mean(data)  # Média dos valores
  total_diff <- sum(outer(data, data, FUN = function(x, y) abs(x - y)))  # Soma das diferenças absolutas
  gini <- total_diff / (2 * n^2 * mean_value)  # Cálculo do Gini
  return(gini)
}

# Calculando o coeficiente de Gini para cada site
gini_results <- data %>%
  group_by(site) %>%                # Agrupa os dados por site
  summarize(gini = gini_coefficient(dap))  # Calcula o Gini por grupo

# Exibindo os resultados
print(gini_results)
gini_results.df <- as.data.frame(gini_results)

# Salvando os resultados

library(writexl)
write_xlsx(gini_results.df,"C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/gini.xlsx")

