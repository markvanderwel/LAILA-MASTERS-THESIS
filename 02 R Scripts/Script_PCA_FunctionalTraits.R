
### PCA FUNCTIONAL TRAITS ###


funcional_regua <- read.csv("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/funcional_regua.csv", 
                            row.names = 1, sep = ";")



funcional_ml <- read.csv("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/funcional_ml.csv", 
                         row.names = 1, sep = ";")

### Incomplete data ###

funcional_regua <- funcional_regua[, -((ncol(funcional_regua)-2):ncol(funcional_regua))]

funcional_ml <- funcional_ml[, -((ncol(funcional_ml)-3):ncol(funcional_ml))]


funcional_regua$specie <- rownames(funcional_regua)
funcional_ml$specie    <- rownames(funcional_ml)

funcional_regua$Site <- "REGUA"
funcional_ml$Site    <- "Mico_Leao"

## Bind the dfs

funcional_total <- rbind(funcional_regua, funcional_ml)

## Run PCA

library(ggplot2)
library(ggfortify)
library(dplyr)

# Select only numerical columns (traits)
traits_num <- funcional_total[, sapply(funcional_total, is.numeric)]

g_pca_simple <- ggplot(scores, aes(x = PC1, y = PC2)) +
  # Eixos tracejados
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
  # Pontos das espécies
  geom_point(size = 2.8, alpha = 0.9, colour = "black") +
  # Vetores dos traits (setas)
  geom_segment(data = loadings,
               aes(x = 0, y = 0, xend = PC1*3, yend = PC2*3),
               arrow = arrow(length = unit(0.25, "cm")),
               colour = "black", linewidth = 1.1) +
  # Rótulos dos traits — LDMC deslocado um pouco mais à direita e acima
  geom_text(data = loadings,
            aes(x = ifelse(Trait == "LDMC", PC1*3.4, PC1*3.2),
                y = ifelse(Trait == "LDMC", PC2*3.4, PC2*3.2),
                label = Trait),
            colour = "black", size = 4, fontface = "bold", hjust = 0) +
  # Eixos proporcionais + margem expandida
  coord_equal(clip = "off") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.margin = margin(10, 35, 10, 10) # folga extra à direita
  ) +
  labs(
    title = "PCA of functional traits",
    subtitle = "Points = species; Arrows = trait loadings",
    x = lab_x,
    y = lab_y
  )

g_pca_simple
ggsave("PCA_functional_traits.jpeg", g_pca_simple,
       width = 16, height = 12, units = "cm", dpi = 600)

