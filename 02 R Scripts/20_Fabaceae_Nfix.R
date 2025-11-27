#### Fabaceae

example <- read_csv("~/01 Masters_LA/00 MASTERS-DATA/01 Datasets/01_raw_data/filogenia_total.csv")
View(example)

fab_species <- example %>% 
  filter(family == "Fabaceae") %>% 
  distinct(species)
fab_species # 34

fab <- read.table("~/01 Masters_LA/00 MASTERS-DATA/01 Datasets/01_raw_data/fabacea_nfix.txt",
                  header = TRUE,
                  sep = ",",
                  stringsAsFactors = FALSE)

fab

library(tidyverse)

fab_summary <- fab %>%
  count(specie.nfix) %>%
  mutate(
    type = ifelse(specie.nfix == 1, "N-fixing", "Non-fixing"),
    pct  = round(n / sum(n) * 100, 1)
  )

fab_summary

library(ggplot2)

g_fix <- ggplot(fab_summary, aes(x = type, y = n)) +
  geom_col(fill = "grey70") +
  geom_text(aes(label = n), vjust = -0.5, size = 5) +
  theme_classic(base_size = 14) +
  labs(
    title = "Number of nitrogen-fixing vs. non-fixing Fabaceae species",
    x = "",
    y = "Number of species"
  )

g_fix

ggsave(
  "~/01 Masters_LA/06 Figures/02 plots/fabaceae_nfix_barplot.jpeg",
  g_fix, width = 6, height = 4, dpi = 300
)
