  - Use built in datasets
  - Usei Tidyverse package to:
    - Select variables
    - Filter observations
    - Create a new variable
    - Create a summary
  
library(tidyverse)
  
View(starwars)

starwars %>% #within
  select(gender,mass,height,species)%>%
  filter(species=="Human")%>%
  na.omit() %>%
  mutate(height=height/100) %>% #quero passar de cm para metros
  mutate(BMI = mass/height^2)%>%
  group_by(gender) %>%
  summarise(Average_BMI= mean(BMI))


# Adicionar a coluna 'nova_coluna' de df2 a df1 sem correspondência por id
meusdados <- meusdados %>%
  mutate(produtividade_gm2ano_c = chave %>% pull(Produtividade_Chave1))
