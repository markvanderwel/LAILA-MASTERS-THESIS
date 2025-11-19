####################################
####### STRUCTURAL DIVERSITY #######
####################################

# Pacote necessário para manipulação de dados
if (!require("dplyr")) install.packages("dplyr")
library(dplyr)
library(readxl)

data <- read_excel("C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/TRABALHO EXCEL/dap.xlsx")

# Função para calcular o coeficiente de Gini
gini_coefficient <- function(data) {
  data <- sort(data)  # Ordenar os dados
  n <- length(data)  # Número de observações
  mean_value <- mean(data)  # Média dos valores
  total_diff <- sum(outer(data, data, FUN = function(x, y) abs(x - y)))  # Soma das diferenças absolutas
  gini <- total_diff / (2 * n^2 * mean_value)  # Cálculo do Gini
  return(gini)
}



# Calculando o coeficiente de Gini para cada site
gini_results <- data %>%
  group_by(site) %>%                # Agrupa os dados por site
  summarize(gini = gini_coefficient(dap))  # Calcula o Gini por grupo

# Exibindo os resultados
print(gini_results)
gini_results.df <- as.data.frame(gini_results)

# Salvando os resultados

library(writexl)
write_xlsx(gini_results.df,"C:/Users/laila/OneDrive/Documentos/2. Mestrado/2. Análise Estatística/gini.xlsx")
