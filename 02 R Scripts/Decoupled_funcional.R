### DECOUPLE ### ML ###

traits <- read.csv("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/funcional_ml.csv", 
                   row.names = 1, sep = ";")

install.packages("cluster")
library(cluster)

Fdist <- daisy(traits, metric = "gower")


# Phylogenetic Tree
example <- read.csv("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/filogenia_ml.csv")

### generate a phylogeny for the sample species list
tree <- phylo.maker(example, tree = GBOTB.extended.TPL,output.sp.list = TRUE,nodes = nodes.info.1.TPL, scenarios="S3")
tree_ok <- tree$scenario.3

Pdist = cophenetic(tree_ok)

# Transform Pdist to eigenvectors trough PCoA

library(vegan)

pcoa_phylo <- cmdscale(Pdist, eig = TRUE, k = nrow(Pdist) - 1)
phylo_axes <- pcoa_phylo$points

# RDA: traits as response, phylogenie as predictor

traits_rda <- rda(traits ~ ., data = as.data.frame(phylo_axes))

# Extract residual scores from RDA (traits decoupled from Phylogeny)

dcF <- scores(traits_rda, display = "sites", choices = 1:traits_rda$CA$rank)
dcFdist <- dist(dcF)  # matriz funcional decouplada da filogenia


## Transform dcFdist into predict variables

n <- attr(dcFdist, "Size")
pcoa_dcF <- cmdscale(dcFdist, eig = TRUE, k = n - 1)
dcF_axes <- as.data.frame(pcoa_dcF$points)

library(writexl)
setwd("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/Resultados funcional")
write_xlsx(dcF_axes,"dcFaxes_ml_reg.xlsx")

pcoa_dcF$eig / sum(pcoa_dcF$eig)  # proporção explicada por cada eixo

library(vegan)

# dcF_axes deve ter as mesmas espécies como rownames que estão nas colunas da matriz_comunidade
comunidade_ml <- read.csv("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/comunidade_ml.csv", 
                          row.names = 1, header = TRUE, sep = ";", check.names = FALSE)

rel_abund <- decostand(comunidade_ml, method = "total")  # abundância relativa
CWM_dcF <- as.matrix(rel_abund) %*% as.matrix(dcF_axes)
CWM_dcF <- as.data.frame(CWM_dcF)
colnames(CWM_dcF) <- paste0("PCoA", seq_len(ncol(CWM_dcF)))  # nomeando os eixos
write_xlsx(CWM_dcF,"CWM_dcF_ml_reg.xlsx")

### DECOUPLE – PHYLOGENY DECOUPLED FROM TRAITS ###

# Phylogeny decoupled from traits

phylo_rda <- rda(phylo_axes ~ ., data = as.data.frame(traits))
dcP <- scores(phylo_rda, display = "sites", choices = 1:traits_rda$CA$rank)
dcPdist <- dist(dcP)

# Load required packages
library(vegan)
library(writexl)

# 1. Perform PCoA on decoupled phylogenetic distance matrix
n_p <- attr(dcPdist, "Size")
pcoa_dcP <- cmdscale(dcPdist, eig = TRUE, k = n_p - 1)
dcP_axes <- as.data.frame(pcoa_dcP$points)

# Standardize column names in the community matrix
colnames(comunidade_ml) <- gsub(" ", "_", colnames(comunidade_ml))

# 2. Check and align species order between community and PCoA axes
# (Assumes 'comunidade_ml' has species as columns and sites as rows)
if (!all(colnames(comunidade_ml) %in% rownames(dcP_axes))) {
  stop("Species names in community matrix and dcP axes do not match.")
}
dcP_axes <- dcP_axes[colnames(comunidade_ml), , drop = FALSE]

setdiff(colnames(comunidade_ml), rownames(dcP_axes))


# 3. Calculate relative abundance
rel_abund <- decostand(comunidade_ml, method = "total")

# 4. Calculate Community-Weighted Means of decoupled phylogenetic axes
CWM_dcP <- as.matrix(rel_abund) %*% as.matrix(dcP_axes)
CWM_dcP <- as.data.frame(CWM_dcP)
colnames(CWM_dcP) <- paste0("PCoA", seq_len(ncol(CWM_dcP)))

# 5. Export to Excel
write_xlsx(CWM_dcP, "CWM_dcP_ml_reg.xlsx")


#
#
#
#
#
#

### DECOUPLE ### REGUA ###

traits <- read.csv("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/funcional_regua.csv", 
                   row.names = 1, sep = ";")

install.packages("cluster")
library(cluster)

Fdist <- daisy(traits, metric = "gower")


# Phylogenetic Tree
example <- read.csv("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/filogenia_regua.csv")

### generate a phylogeny for the sample species list
tree <- phylo.maker(example, tree = GBOTB.extended.TPL,output.sp.list = TRUE,nodes = nodes.info.1.TPL, scenarios="S3")
tree_ok <- tree$scenario.3

Pdist = cophenetic(tree_ok)

# Transform Pdist to eigenvectors trough PCoA

library(vegan)

pcoa_phylo <- cmdscale(Pdist, eig = TRUE, k = nrow(Pdist) - 1)
phylo_axes <- pcoa_phylo$points

# RDA: traits as response, phylogenie as predictor

traits_rda <- rda(traits ~ ., data = as.data.frame(phylo_axes))

# Extract residual scores from RDA (traits decoupled from Phylogeny)

dcF <- scores(traits_rda, display = "sites", choices = 1:traits_rda$CA$rank)
dcFdist <- dist(dcF)  # matriz funcional decouplada da filogenia

## Transform dcFdist into predict variables

n <- attr(dcFdist, "Size")
pcoa_dcF <- cmdscale(dcFdist, eig = TRUE, k = n - 1)
dcF_axes <- as.data.frame(pcoa_dcF$points)

library(writexl)
setwd("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/Resultados funcional")

library(vegan)

comunidade_regua <- read.csv("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/comunidade_regua.csv", 
                          row.names = 1, header = TRUE, sep = ";", check.names = FALSE)

rel_abund <- decostand(comunidade_regua, method = "total")  # abundância relativa
CWM_dcF <- as.matrix(rel_abund) %*% as.matrix(dcF_axes)
CWM_dcF <- as.data.frame(CWM_dcF)
colnames(CWM_dcF) <- paste0("PCoA", seq_len(ncol(CWM_dcF)))  # nomeando os eixos
write_xlsx(CWM_dcF,"CWM_dcF_regua_reg.xlsx")

### DECOUPLE – PHYLOGENY DECOUPLED FROM TRAITS ###

# Phylogeny decoupled from traits

phylo_rda <- rda(phylo_axes ~ ., data = as.data.frame(traits))
dcP <- scores(phylo_rda, display = "sites", choices = 1:traits_rda$CA$rank)
dcPdist <- dist(dcP)

# Load required packages
library(vegan)
library(writexl)

# 1. Perform PCoA on decoupled phylogenetic distance matrix
n_p <- attr(dcPdist, "Size")
pcoa_dcP <- cmdscale(dcPdist, eig = TRUE, k = n_p - 1)
dcP_axes <- as.data.frame(pcoa_dcP$points)

# Standardize column names in the community matrix
colnames(comunidade_regua) <- gsub(" ", "_", colnames(comunidade_regua))

# 2. Check and align species order between community and PCoA axes
# (Assumes 'comunidade_ml' has species as columns and sites as rows)
if (!all(colnames(comunidade_regua) %in% rownames(dcP_axes))) {
  stop("Species names in community matrix and dcP axes do not match.")
}
dcP_axes <- dcP_axes[colnames(comunidade_regua), , drop = FALSE]

setdiff(colnames(comunidade_regua), rownames(dcP_axes))


# 3. Calculate relative abundance
rel_abund <- decostand(comunidade_regua, method = "total")

# 4. Calculate Community-Weighted Means of decoupled phylogenetic axes
CWM_dcP <- as.matrix(rel_abund) %*% as.matrix(dcP_axes)
CWM_dcP <- as.data.frame(CWM_dcP)
colnames(CWM_dcP) <- paste0("PCoA", seq_len(ncol(CWM_dcP)))

# 5. Export to Excel
write_xlsx(CWM_dcP, "CWM_dcP_regua_reg.xlsx")
