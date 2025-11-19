### DADOS AMBIENTAIS ###

# BAIXAR NO TERRA CLIMATE

#-/-

# Raster
# obs: ppt = precipitation (mm); pet = potential_evaporation (mm)
## tmax (média mensal) = temperatura máxima (ºC); tmin (média mensal) = temperatura mínima (ºC)
### vpd (média mensal) = Vapor Pressure Deficit (kpd)

install.packages("ncdf4")
library(ncdf4)
library(raster)

##Carrega o raster como stack
# ppt
terraclimate_rst <- stack(list.files("C:/Users/laila/OneDrive/Documentos/2. Mestrado/3. R/dados_ambientais", pattern = "TerraClimate_ppt", full.names = T))
#tmax
terraclimate_rst <- stack(list.files("C:/Users/laila/OneDrive/Documentos/2. Mestrado/3. R/dados_ambientais", pattern = "TerraClimate_tmax", full.names = T))
#tmin
terraclimate_rst <- stack(list.files("C:/Users/laila/OneDrive/Documentos/2. Mestrado/3. R/dados_ambientais", pattern = "TerraClimate_tmin", full.names = T))
# pet
terraclimate_rst <- stack(list.files("C:/Users/laila/OneDrive/Documentos/2. Mestrado/3. R/dados_ambientais", pattern = "TerraClimate_pet", full.names = T))
# vpd
terraclimate_rst <- stack(list.files("C:/Users/laila/OneDrive/Documentos/2. Mestrado/3. R/dados_ambientais", pattern = "TerraClimate_vpd", full.names = T))
warnings()

##renomear os layers do stack
names(terraclimate_rst) <- apply(expand.grid(1:12, 2018:2022),1, FUN = paste, collapse = "_")

##data.frame com as coordenadas latlong

library(readxl)
moderadores <- read_excel("C:/Users/laila/OneDrive/Documentos/2. Mestrado/Mapa/moderadores.xlsx")
View(moderadores)

##Cortar o raster em um extent para reduzir o uso da memória
e1 <- extent(-44.608007, -40.974161, -22.859079, -21.272083)
terraclimate_rst_crop <- crop(terraclimate_rst, e1)

##extrair os dados climáticos
env_data <- extract(terraclimate_rst_crop, moderadores[, c("longitude", "latitude")])

## Salvar em planilha do excel

install.packages("writexl")
library(writexl)
env_data_df <- as.data.frame(env_data)
write_xlsx(env_data_df, "C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Dados/amb_vpd.xlsx")
