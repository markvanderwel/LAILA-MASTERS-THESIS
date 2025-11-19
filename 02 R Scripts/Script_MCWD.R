library(readxl) # Para carregar arquivos CSV
library(sf)    # Para lidar com dados espaciais

# Carregar a planilha
sites <- read_excel("C:/Users/laila/OneDrive/Documentos/2. Mestrado/6. Mapa/moderadores.xlsx")

# Converter para objeto espacial
sites_sf <- st_as_sf(sites, coords = c("Longitude", "Latitude"), crs = 4326) # CRS WGS84

library(terra)

### testar se está tudo certo

# Carregar o raster de 2016 (como exemplo)
mcwd_2016 <- rast("C:/Global_MCWD/MCWD_2016.tif")

# Verificar a extensão
print(ext(mcwd_2016))

# Visualizar o raster
plot(mcwd_2016, main = "MCWD 2016")

# Extrair valores para os pontos
mcwd_values_2016 <- extract(mcwd_2016, vect(sites_sf), method = "nearest")

# Verificar o resultado
print(head(mcwd_values_2016))

# Plotar o raster
plot(mcwd_2016, main = "MCWD 2016")

# Adicionar os pontos
plot(st_geometry(sites_sf), add = TRUE, col = "red", pch = 20, cex = 1.5)

#############################################################################
# TODOS OS SITES DE UMA VEZ

# Listar os arquivos MCWD de 2016 a 2020
raster_files <- list.files(
  path = "C:/Global_MCWD",
  pattern = "MCWD_201[6-9]|MCWD_2020\\.tif$", 
  full.names = TRUE
)

# Carregar todos os arquivos como um 'SpatRaster'
mcwd_stack <- rast(raster_files)

names(mcwd_stack) <- c("MCWD_2016", "MCWD_2017", "MCWD_2018", "MCWD_2019", "MCWD_2020")
print(mcwd_stack) # conferindo se está certo

sites_sf <- st_transform(sites_sf, crs(mcwd_stack))

# Extrair valores para todos os sites e anos
mcwd_values <- extract(mcwd_stack, vect(sites_sf))

# Renomear as colunas para identificar os anos
colnames(mcwd_values) <- c("ID", paste0("MCWD_", 2016:2020))

# Combinar os valores extraídos com a tabela original
sites <- cbind(sites, mcwd_values[, -1]) # Excluir a coluna ID

install.packages("writexl")

library(writexl)

write_xlsx(sites, "C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/mcwd_sites.xlsx")




















