## ======================================================
##  Summary statistics for the first Results paragraph
##  Mean ± SE, Min–Max, and PCPS variance (%)
## ======================================================

library(writexl)
library(readxl)
setwd("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/TRABALHO EXCEL")
dadosmisto <- read_excel("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/dadosmisto.xlsx")

# funções auxiliares
se_ <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
summ_ <- function(x) c(
  n = sum(!is.na(x)),
  mean = mean(x, na.rm = TRUE),
  se = se_(x),
  min = min(x, na.rm = TRUE),
  max = max(x, na.rm = TRUE)
)
fmt_pm <- function(mean, se, digits = 2) paste0(round(mean, digits), " ± ", round(se, digits))
fmt_range <- function(minv, maxv, digits = 2) paste0(round(minv, digits), "–", round(maxv, digits))

# variáveis de interesse
vars <- list(
  `Net carbon assimilation (g m^-2 yr^-1)` = dadosmisto1$produt_z_g.ano,
  `Species richness (SR)`                  = dadosmisto1$SR,
  `Phylogenetic diversity (PD)`            = dadosmisto1$SESPDab,
  `Phylogenetic composition (PCPS1)`       = dadosmisto1$pcps1ab,
  `Silt (%)`                               = dadosmisto1$silte
)

# calcular estatísticas
res_list <- lapply(vars, summ_)
res_df <- as.data.frame(do.call(rbind, res_list))
res_df$`Mean ± SE` <- mapply(fmt_pm, res_df$mean, res_df$se)
res_df$Range <- mapply(fmt_range, res_df$min, res_df$max)
res_out <- res_df[, c("n","Mean ± SE","Range","mean","se","min","max")]
rownames(res_out) <- names(vars)

cat("\n==== Summary for Results paragraph ====\n")
print(res_out)

res_out_export <- res_out
res_out_export <- cbind(
  Metric = rownames(res_out_export),
  res_out_export
)
rownames(res_out_export) <- NULL

write_xlsx(as.data.frame(res_out),
           "C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/TRABALHO EXCEL/summary_results_paragraph_withCSA.xlsx")


