ntt_pcps <- pcps(ntt_phylodata$community[,] , ntt_phylodata$phylodist)

scores_pcps <- scores.pcps(ntt_pcps, choice = c(1,5))

scores_pcps_plots <- scores_pcps$scores.sites[rownames(scores_pcps$scores.sites) %in% data_sem$X,]
scores_pcps_spp <- scores_pcps$scores.species

apg_arrows <- apply(scores_pcps_spp,2, function(x) tapply(x, list(apg_groups), mean)) #apg_groups é o vetor das famílias de cada espécie

library(ggplot2)

color_scale <- scales::gradient_n_pal(c("orange", "darkgreen"))(range(data_sem$deciduousness))

apg_arrows <- as.data.frame(apg_arrows)

apg_arrows$angle <- atan2(apg_arrows$pcps.5, apg_arrows$pcps.1) * (180/pi)

apg_arrows$labels <- c("Anac", "Euph", "Faba", "Laur", "Mela", "Meli", "Mora", "Myrt", "Rubi", "Sali")


pdf("fig_pcps.pdf", height = 6, width = 7)

ggplot(as.data.frame(scores_pcps_plots), aes(x = pcps.1, y = pcps.5)) +
  geom_point(aes(color = data_sem$deciduousness)) +
  scale_color_gradientn(colours = color_scale) +
  labs(color = "Deciduousness", x = "PCPS1", y = "PCPS5") +
  geom_text(data = as.data.frame(apg_arrows),
            aes(x = pcps.1, y = pcps.5, label = labels, fontface = "bold"), color = "blue", size = 5) +
  geom_label(data = as.data.frame(apg_arrows),
             aes(x = pcps.1, y = pcps.5, label = labels), 
             fontface = "bold", color = "black", size = 5,
             fill = "grey", alpha = 0.5) +
  theme_minimal() + 
  guides(color = FALSE) +
  theme(axis.line = element_line(color = "black", size = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))
dev.off()
