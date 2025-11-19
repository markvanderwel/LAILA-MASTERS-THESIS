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

