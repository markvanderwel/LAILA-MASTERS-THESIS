## EXPORTANDO DADOS DO TRY ##

install.packages("rtry")
library(rtry)

# Importar o arquivo
input <- rtry_import("C:/Users/laila/OneDrive/Documentos/2. Mestrado/TRY/33506.txt")

#Criar um dataframe apenas com as colunas de interesse
data1 <- rtry_select_col(input,TraitID,SpeciesName,DataID,StdValue,OrigObsDataID,UnitName,ObservationID,Reference)



#Selecionar as linhas que contém a densidade da madeira (4)
#'Transformar em dataframe
dataf <- as.data.frame(data1)
#' Load dplyr package
library(dplyr)
#' Using dplyr::filter
DM <- dplyr::filter(dataf, TraitID %in% c("4"))

#Remover colunas que não me interessam
DM1 <- rtry_remove_col(DM,OrigObsDataID,ObservationID)

setwd("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Dados")
write.csv(DM1, file = "DM_Try")


