teste <- read.csv("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Dados/TRABALHO EXCEL/especie.regua.csv",header = T, sep = ";")

especie_total <- teste$especie_total
especie_dom <- teste$especie_dom

especies_comuns <- intersect(especie_total, especie_dom)

regua_comum <- as.data.frame(especies_comuns)

setwd("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Dados/TRABALHO EXCEL")
write.csv(regua_comum, file = "REGUA_comum.csv")