
####################
### TRY FABACEAE ###
####################

install.packages("rtry")
library(rtry)
library(dplyr)
library(stringr)

try_data <- rtry_import("~/01 Masters_LA/00 MASTERS-DATA/01 Datasets/01_raw_data/Try_Fabaceae/45484.txt")

head(try_data)
str(try_data)

unique(try_data$TraitID) # NA is SLA?
unique(try_data$TraitName)

##########
## LDMC ##
##########

# 1) Filter LDMC (TraitID = 47)
ldmc_raw <- try_data %>%
  filter(TraitID == 47) %>%               # LDMC
  select(SpeciesName, StdValue, UnitName)

head(ldmc_raw)

# 2) Remove NAs
ldmc_raw <- ldmc_raw %>%
  filter(!is.na(StdValue))

# 3) Calculate mean per specie
ldmc_species <- ldmc_raw %>%
  group_by(SpeciesName, UnitName) %>%     
  summarise(
    LDMC_mean = mean(StdValue, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  )

head(ldmc_species)

# 3.2) Clean data

ldmc_species_clean <- ldmc_species %>%
  mutate(
    # remove duplicate spaces
    Species_trim = str_squish(SpeciesName),
    # counts the number of words
    n_words = str_count(Species_trim, "\\S+")
  ) %>%
  # It only keeps names with exactly 2 words (Genus + epithet).
  filter(n_words == 2) %>%
  # now standardize the final name
  mutate(Species_std = Species_trim)

library(ggplot2)

ldmc_plot <- ggplot(ldmc_species_clean, aes(x = "", y = LDMC_mean)) +
  geom_boxplot(fill = "grey80", color = "black", width = 0.4) +
  geom_jitter(width = 0.1, alpha = 0.5, color = "darkorange") +
  theme_classic(base_size = 14) +
  labs(
    title = "Distribution of LDMC (TRY Trait 47)",
    y = "LDMC (mean StdValue)",
    x = ""
  )

ldmc_plot

########
## WD ## Good to compare
########

# 1) Filter WD (TraitID = 4)
wd_raw <- try_data %>%
  filter(TraitID == 4) %>%               # WD
  select(SpeciesName, StdValue, UnitName)

head(wd_raw)

# 2) Remove NAs
wd_raw <- wd_raw %>%
  filter(!is.na(StdValue))


# 3.1) Calculate mean per specie
wd_species <- wd_raw %>%
  group_by(SpeciesName, UnitName) %>%     
  summarise(
    wd_mean = mean(StdValue, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  )

head(wd_species)

# 3.2) Clean data

wd_species_clean <- wd_species %>%
  mutate(
    # remove duplicate spaces
    Species_trim = str_squish(SpeciesName),
    # counts the number of words
    n_words = str_count(Species_trim, "\\S+")
  ) %>%
  # It only keeps names with exactly 2 words (Genus + epithet).
  filter(n_words == 2) %>%
  # now standardize the final name
  mutate(Species_std = Species_trim)


library(ggplot2)

wd_plot <-  ggplot(wd_species_clean, aes(x = "", y = wd_mean)) +
  geom_boxplot(fill = "grey80", color = "black", width = 0.4) +
  geom_jitter(width = 0.1, alpha = 0.5, color = "darkorange") +
  theme_classic(base_size = 14) +
  labs(
    title = "Distribution of WD (TRY Trait 4)",
    y = "WD (mean StdValue)",
    x = ""
  )

wd_plot

####################
## Leaf N content ##
####################

# 1) Filter LNC (TraitID = 50)
lnc_raw <- try_data %>%
  filter(TraitID == 50) %>%               # LNC
  select(SpeciesName, StdValue, UnitName)

head(lnc_raw)

# 2) Remove NAs
lnc_raw <- lnc_raw %>%
  filter(!is.na(StdValue))


# 3.1) Calculate mean per specie
lnc_species <- lnc_raw %>%
  group_by(SpeciesName, UnitName) %>%     
  summarise(
    lnc_mean = mean(StdValue, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  )

head(lnc_species)

# 3.2) Clean data

lnc_species_clean <- lnc_species %>%
  mutate(
    # remove duplicate spaces
    Species_trim = str_squish(SpeciesName),
    # counts the number of words
    n_words = str_count(Species_trim, "\\S+")
  ) %>%
  # It only keeps names with exactly 2 words (Genus + epithet).
  filter(n_words == 2) %>%
  # now standardize the final name
  mutate(Species_std = Species_trim)


library(ggplot2)

lnc_plot <- ggplot(lnc_species_clean, aes(x = "", y = lnc_mean)) +
  geom_boxplot(fill = "grey80", color = "black", width = 0.4) +
  geom_jitter(width = 0.1, alpha = 0.5, color = "darkorange") +
  theme_classic(base_size = 14) +
  labs(
    title = "Distribution of Leaf N content (TRY Trait 50)",
    y = "LNC (mean StdValue)",
    x = ""
  )

lnc_plot

####################
## Leaf P content ##
####################

# 1) Filter LNC (TraitID = 50)
lpc_raw <- try_data %>%
  filter(TraitID == 51) %>%               # LPC
  select(SpeciesName, StdValue, UnitName)

head(lpc_raw)

# 2) Remove NAs
lpc_raw <- lpc_raw %>%
  filter(!is.na(StdValue))


# 3.1) Calculate mean per specie
lpc_species <- lpc_raw %>%
  group_by(SpeciesName, UnitName) %>%     
  summarise(
    lpc_mean = mean(StdValue, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  )

head(lpc_species)

# 3.2) Clean data

lpc_species_clean <- lpc_species %>%
  mutate(
    # remove duplicate spaces
    Species_trim = str_squish(SpeciesName),
    # counts the number of words
    n_words = str_count(Species_trim, "\\S+")
  ) %>%
  # It only keeps names with exactly 2 words (Genus + epithet).
  filter(n_words == 2) %>%
  # now standardize the final name
  mutate(Species_std = Species_trim)


library(ggplot2)

lpc_plot <- ggplot(lpc_species_clean, aes(x = "", y = lpc_mean)) +
  geom_boxplot(fill = "grey80", color = "black", width = 0.4) +
  geom_jitter(width = 0.1, alpha = 0.5, color = "darkorange") +
  theme_classic(base_size = 14) +
  labs(
    title = "Distribution of Leaf P content (TRY Trait 51)",
    y = "LPC (mean StdValue)",
    x = ""
  )

lpc_plot

###########
### SLA ###
###########

try_data_SLA <- rtry_import("~/01 Masters_LA/00 MASTERS-DATA/01 Datasets/01_raw_data/Try_Fabaceae/sla.txt")

unique(try_data_SLA$TraitID)

sla_count <- try_data_SLA %>%
  group_by(TraitID) %>%
  summarise(n_obs = n())

# ID 3116 is the chosen

# 1) Filter SLA (TraitID = 3116)
sla_raw <- try_data_SLA %>%
  filter(TraitID == 3116) %>%               # SLA
  select(SpeciesName, StdValue, UnitName)

head(sla_raw)

# 2) Remove NAs
sla_raw <- sla_raw %>%
  filter(!is.na(StdValue))


# 3.1) Calculate mean per specie
sla_species <- sla_raw %>%
  group_by(SpeciesName, UnitName) %>%     
  summarise(
    sla_mean = mean(StdValue, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  )

head(sla_raw)

# 3.2) Clean data

sla_species_clean <- sla_species %>%
  mutate(
    # remove duplicate spaces
    Species_trim = str_squish(SpeciesName),
    # counts the number of words
    n_words = str_count(Species_trim, "\\S+")
  ) %>%
  # It only keeps names with exactly 2 words (Genus + epithet).
  filter(n_words == 2) %>%
  # now standardize the final name
  mutate(Species_std = Species_trim)


library(ggplot2)

sla_plot <- ggplot(sla_species_clean, aes(x = "", y = sla_mean)) +
  geom_boxplot(fill = "grey80", color = "black", width = 0.4) +
  geom_jitter(width = 0.1, alpha = 0.5, color = "darkorange") +
  theme_classic(base_size = 14) +
  labs(
    title = "Distribution of SLA (TRY Trait 3116)",
    y = "SLA (mean StdValue)",
    x = ""
  )

sla_plot

#########################
####### OTHER SPP #######
#########################

try_data_other <- rtry_import("~/01 Masters_LA/00 MASTERS-DATA/01 Datasets/01_raw_data/Try_Fabaceae/spp_pos.txt")

unique(try_data_other$TraitID)

##########
## LDMC ##
##########

# 1) Filter LDMC (TraitID = 47)
ldmc_raw1 <- try_data_other %>%
  filter(TraitID == 47) %>%               # LDMC
  select(SpeciesName, StdValue, UnitName)

head(ldmc_raw)

# 2) Remove NAs
ldmc_raw1 <- ldmc_raw1 %>%
  filter(!is.na(StdValue))

# 3) Calculate mean per specie
ldmc_species1 <- ldmc_raw1 %>%
  group_by(SpeciesName, UnitName) %>%     
  summarise(
    LDMC_mean = mean(StdValue, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  )

head(ldmc_species1)

# 3.2) Clean data

ldmc_species_clean1 <- ldmc_species1 %>%
  mutate(
    # remove duplicate spaces
    Species_trim = str_squish(SpeciesName),
    # counts the number of words
    n_words = str_count(Species_trim, "\\S+")
  ) %>%
  # It only keeps names with exactly 2 words (Genus + epithet).
  filter(n_words == 2) %>%
  # now standardize the final name
  mutate(Species_std = Species_trim)

library(ggplot2)

ldmc_plot1 <- ggplot(ldmc_species_clean1, aes(x = "", y = LDMC_mean)) +
  geom_boxplot(fill = "grey80", color = "black", width = 0.4) +
  geom_jitter(width = 0.1, alpha = 0.5, color = "darkorange") +
  theme_classic(base_size = 14) +
  labs(
    title = "Distribution of LDMC (TRY Trait 47)",
    y = "LDMC (mean StdValue)",
    x = ""
  )

ldmc_plot1

# Join 2 tables to do boxplot

ldmc_fab <- ldmc_species_clean %>%
  mutate(group = "Fabaceae") %>%
  select(Species_std, LDMC_mean, UnitName, group)

ldmc_pos <- ldmc_species_clean1 %>%
  mutate(group = "Other PCPS1 families") %>%
  select(Species_std, LDMC_mean, UnitName, group)

ldmc_all <- bind_rows(ldmc_fab, ldmc_pos)

library(ggplot2)

gg_ldmc_comp <- ggplot(ldmc_all, aes(x = group, y = LDMC_mean)) +
  geom_boxplot(fill = "grey90", color = "black") +
  geom_jitter(width = 0.1, alpha = 0.6, aes(color = group)) +
  theme_classic(base_size = 14) +
  labs(
    title = "LDMC of Fabaceae vs. other PCPS1 families (TRY Trait 47)",
    x = "",
    y = "LDMC (mean StdValue)"
  ) +
  theme(legend.position = "none")

gg_ldmc_comp

t_test_ldmc <- t.test(LDMC_mean ~ group, data = ldmc_all)
t_test_ldmc # p-value = 0.8816

########
## WD ## other families
########

wd_raw1 <- try_data_other %>%
  filter(TraitID == 4) %>%
  select(SpeciesName, StdValue, UnitName)

wd_raw1 <- wd_raw1 %>%
  filter(!is.na(StdValue))

wd_species1 <- wd_raw1 %>%
  group_by(SpeciesName, UnitName) %>%
  summarise(
    wd_mean = mean(StdValue, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  )

wd_species_clean1 <- wd_species1 %>%
  mutate(
    Species_trim = str_squish(SpeciesName),
    n_words = str_count(Species_trim, "\\S+")
  ) %>%
  filter(n_words == 2) %>%
  mutate(Species_std = Species_trim)

# Join

wd_fab <- wd_species_clean %>%
  mutate(group = "Fabaceae") %>%
  select(Species_std, wd_mean, group)

wd_pos <- wd_species_clean1 %>%
  mutate(group = "Other PCPS1 families") %>%
  select(Species_std, wd_mean, group)

wd_all <- bind_rows(wd_fab, wd_pos)

# Plot

gg_wd_comp <- ggplot(wd_all, aes(x = group, y = wd_mean)) +
  geom_boxplot(fill = "grey90", color = "black") +
  geom_jitter(width = 0.1, alpha = 0.6, aes(color = group)) +
  theme_classic(base_size = 14) +
  labs(
    title = "WD of Fabaceae vs. other PCPS1 families (TRY Trait 4)",
    x = "",
    y = "Wood density (mean StdValue)"
  ) +
  theme(legend.position = "none")

gg_wd_comp

t_test_wd <- t.test(wd_mean ~ group, data = wd_all)
t_test_wd # p-value = 0.09385

#########
## SLA ##
#########

# 1) Filter SLA (TraitID = 3116)
sla_raw1 <- try_data_other %>%
  filter(TraitID == 3116) %>%               # SLA
  select(SpeciesName, StdValue, UnitName)

head(sla_raw1)

# 2) Remove NAs
sla_raw1 <- sla_raw1 %>%
  filter(!is.na(StdValue))

# 3) Calculate mean per specie
sla_species1 <- sla_raw1 %>%
  group_by(SpeciesName, UnitName) %>%     
  summarise(
    sla_mean = mean(StdValue, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  )

head(sla_species1)

# 3.2) Clean data

sla_species_clean1 <- sla_species1 %>%
  mutate(
    # remove duplicate spaces
    Species_trim = str_squish(SpeciesName),
    # counts the number of words
    n_words = str_count(Species_trim, "\\S+")
  ) %>%
  # It only keeps names with exactly 2 words (Genus + epithet).
  filter(n_words == 2) %>%
  # now standardize the final name
  mutate(Species_std = Species_trim)

library(ggplot2)

sla_plot1 <- ggplot(sla_species_clean1, aes(x = "", y = sla_mean)) +
  geom_boxplot(fill = "grey80", color = "black", width = 0.4) +
  geom_jitter(width = 0.1, alpha = 0.5, color = "darkorange") +
  theme_classic(base_size = 14) +
  labs(
    title = "Distribution of SLA (TRY Trait 3116)",
    y = "SLA (mean StdValue)",
    x = ""
  )

sla_plot1

# Join

sla_fab <- sla_species_clean %>%
  mutate(group = "Fabaceae") %>%
  select(Species_std, sla_mean, group)

sla_pos <- sla_species_clean1 %>%
  mutate(group = "Other PCPS1 families") %>%
  select(Species_std, sla_mean, group)

sla_all <- bind_rows(sla_fab, sla_pos)

# Plot

gg_sla_comp <- ggplot(sla_all, aes(x = group, y = sla_mean)) +
  geom_boxplot(fill = "grey90", color = "black") +
  geom_jitter(width = 0.1, alpha = 0.6, aes(color = group)) +
  theme_classic(base_size = 14) +
  labs(
    title = "SLA of Fabaceae vs. other PCPS1 families (TRY Trait 4)",
    x = "",
    y = "Specific Leaf Area (mean StdValue)"
  ) +
  theme(legend.position = "none")

gg_sla_comp

t_test_wd <- t.test(wd_mean ~ group, data = wd_all)
t_test_wd # p-value = 0.09385

var.test(sla_mean ~ group, data = sla_all) # p-value = 0.5872

# The effect of PCPS1 on productivity is not simply explained by Fabaceae having lower LDMC, lower WD, or higher SLA than other families.

# This opens up the possibility that other lineage-related mechanisms, especially N fixation, are behind the productivity pattern.



