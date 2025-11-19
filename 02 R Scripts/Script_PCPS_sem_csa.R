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
