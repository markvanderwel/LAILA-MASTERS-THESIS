library(ggplot2)
library(readxl)
library(tidyverse)

# Ler os dados
granulometria <- read_excel("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/granulometria.xlsx")

# Remover a coluna 'site' para fazer o PCA com apenas as variáveis numéricas
dados_granulometria <- granulometria[, -1]

# Realizar o PCA usando a função base prcomp()
pca_resultado <- prcomp(dados_granulometria, scale. = TRUE)

# Visualizar o resumo do PCA (variância explicada por cada componente)
summary(pca_resultado)

# Plotar o Screeplot para ver a variância explicada por cada componente
screeplot(pca_resultado, main = "Screeplot - PCA")

# Obter os escores dos componentes principais (PCs)
escores_pca <- pca_resultado$x

# Visualizar as primeiras linhas dos escores
head(escores_pca)

# Adicionar os escores PCA ao dataframe original
granulometria_pca <- cbind(granulometria, escores_pca)

# Selecionar os eixos mais explicativos (exemplo: PC1 e PC2)
granulometria_pca <- granulometria_pca %>%
  select(site, PC1, PC2)

# Exemplo de regressão com os eixos principais mais explicativos
# Ajustando o modelo linear com os componentes principais
modelo_pca <- lm(resposta ~ PC1 + PC2, data = granulometria_pca)

# Resumo do modelo de regressão
summary(modelo_pca)

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

