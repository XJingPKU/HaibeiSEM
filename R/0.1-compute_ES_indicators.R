# Generate indicators of ecosystem services
# XJ
# last edited in 2024/04/12

# Objectives
# 1. Calculate forage production
# 2. Calculate forage quality (community weighted mean of leaf nitrogen)
# 3. Calculate species richness and beta-diversity (LCBD)
# 4. Calculate soil carbon storage
# 5. Calculate soil respiration

###########################################################
# Load library ------------------------
library(tidyverse)
library(FD)
library(adespatial)

rm(list = ls())
# Load data ---------------------------
# experiment design
exp_design <- read.csv("./data/experiment_design.csv")
# species table
spp_tab <- read.csv("./data/plant_species_table-2022-05-28.csv")
spp_code <- read.csv("./data/plant_species_code-2022-05-28.csv")
# leaf nitrogen
leaf_N_early <- read.csv("./data/2012-2013_leaf_nitrogen.csv")
leaf_N_middle <- read.csv("./data/2017_leaf_nitrogen.csv")
leaf_N_late <- read.csv("./data/2020_leaf_nitrogen.csv")
# soil total carbon
soil_C_early <- read.csv("./data/2013_soil_total_carbon.csv")
soil_C_middle <- read.csv("./data/2015_soil_total_carbon.csv")
soil_C_late <- read.csv("./data/2020_soil_total_carbon.csv")
soil_BD <- read.csv("./data/2020_bulk_density.csv")
# soil respiration
year <- c(2012, 2013, 2015, 2016, 2017, 2020)
plts <- paste("plot", c(1:6, 31:36), sep = "")

sr_8100 <- sr_8150 <- list()

for (i in 1:length(year)) {
  path_8100 <- glue::glue("./data/soil_respiration/{year[i]}_soil_respiration_8100.csv")
sr_8100[[i]] <- read.csv(glue::glue(path_8100))
}


for (i in 1:length(plts)) {
  path_8150 <- glue::glue("./data/soil_respiration/{plts[i]}_soil_respiration_8150.csv")
  sr_8150[[i]] <- read.csv(glue::glue(path_8150))
}


###########################################################
# Clean data --------------------------
# species table
spp_tab_clean <- spp_tab %>% 
  rename(idd = N) %>% 
  pivot_longer(cols = Lonicera_rupicola:Carex_pygmaea,
               names_to = "species_name",
               values_to = "biomass") %>% 
  group_by(idd) %>% 
  filter(sum(biomass) > 0) %>% 
  ungroup() %>% 
  group_by(species_name) %>% 
  filter(sum(biomass) > 0) %>% 
  ungroup() 

# leaf nitrogen
leaf_N_early <- leaf_N_early %>% 
  select(species, LNC)
leaf_N_middle <- leaf_N_middle %>% 
  select(species, LNC)
leaf_N_late <- leaf_N_late %>% 
  select(species, LNC)

leaf_N <- leaf_N_early %>% 
  bind_rows(leaf_N_middle) %>% 
  bind_rows(leaf_N_late) %>% 
  na.omit() %>% 
  group_by(species) %>% 
  summarise(leaf_N = mean(LNC)) %>% 
  mutate(species = ifelse(species == "Medicago_archiducis-nicolaii",
                          "Medicago_archiducis.nicolaii", species)) %>% 
  as.data.frame()

# soil total carbon
soil_C_early <- soil_C_early %>% 
  rename(block = Block,
         plot = Plot,
         layer = Layer,
         STC = STC_2013_n) %>% 
  select(block, plot, layer, STC) %>% 
  mutate(exp_period = "Early",
         layer = factor(layer,
                        levels = c(1:5),
                        labels = c("0-10cm", "10-20cm", "20-30cm",
                                   "30-40cm", "40-70cm")))
# replace the missing value of STC with that measured in 2015
soil_C_early$STC[which(is.na(soil_C_early$STC))] <- 84.33


soil_C_middle <- soil_C_middle %>% 
  rename(block = Block,
         plot = Plot,
         layer = Layer,
         STC = STC_2015_n) %>% 
  select(block, plot, layer, STC) %>% 
  mutate(exp_period = "Middle",
         layer = factor(layer,
                        levels = c(1:5),
                        labels = c("0-10cm", "10-20cm", "20-30cm",
                                   "30-40cm", "40-70cm")))

soil_C_late <- soil_C_late %>% 
  select(block, plot, 
         X0.10cm, X10.20cm, X20.30cm, X30.40cm,
         X40.50cm, X50.60cm, X60.70cm) %>% 
  pivot_longer(cols = X0.10cm:X60.70cm,
               names_to = "layer",
               values_to = "STC") %>% 
  mutate(exp_period = "Late",
         layer = factor(layer,
                        levels = c("X0.10cm", "X10.20cm", "X20.30cm", 
                                   "X30.40cm", "X40.50cm", "X50.60cm", 
                                   "X60.70cm"),
                        labels = c("0-10cm", "10-20cm", "20-30cm",
                                   "30-40cm", "40-50cm", "50-60cm",
                                   "60-70cm")))
soil_BD <- soil_BD %>% 
  select(block, plot, 
         X0.10cm, X10.20cm, X20.30cm, X30.40cm,
         X40.50cm, X50.60cm, X60.70cm) %>% 
  pivot_longer(cols = X0.10cm:X60.70cm,
               names_to = "layer",
               values_to = "BD") %>% 
  mutate(layer = factor(layer,
                        levels = c("X0.10cm", "X10.20cm", "X20.30cm", 
                                   "X30.40cm", "X40.50cm", "X50.60cm", 
                                   "X60.70cm"),
                        labels = c("0-10cm", "10-20cm", "20-30cm",
                                   "30-40cm", "40-50cm", "50-60cm",
                                   "60-70cm")))

# soil respiration

for (i in 1:length(year)) {
  if (i %in% c(1:5)) {
    if (i == 4) {
      sr_8100[[i]] <- sr_8100[[i]][1:12]
      names(sr_8100[[i]]) <- tolower(names(sr_8100[[i]]))
      sr_8100[[i]] <- sr_8100[[i]] %>%
        rename(csr = 'csr.n') %>% 
        select(block, plot, year, month, day, csr) %>%
        filter(month %in% c(6:9),
               block %in% c("C", "D", "E"),
               !plot %in% c("W1", "W2", "W3"))
    } else {
      names(sr_8100[[i]]) <- tolower(names(sr_8100[[i]]))
      sr_8100[[i]] <- sr_8100[[i]] %>% 
        select(block, plot, year, month, day, csr) %>% 
        filter(month %in% c(6:9),
               block %in% c("C", "D", "E"),
               !plot %in% c("W1", "W2", "W3"),
               !plot %in% c("w1", "w2", "w3"))
    }
  } else {
    names(sr_8100[[i]]) <- tolower(names(sr_8100[[i]]))
    sr_8100[[i]] <- sr_8100[[i]] %>% 
      mutate(date = as.character(date),
             month = format(as.Date(date, format = "%m/%d/%Y"), "%m"),
             day = format(as.Date(date, format = "%m/%d/%Y"), "%d")) %>% 
      mutate(month = as.numeric(month),
             day = as.numeric(day)) %>% 
      rename(csr = sr) %>% 
      select(block, plot, year, month, day, csr) %>% 
      filter(month %in% c(6:9))
  }
}
sr_8100 <- do.call(rbind, sr_8100) %>% 
  mutate(date = paste(year, month, day, sep = "-")) %>% 
  as.data.frame() %>% 
  mutate(plot = as.numeric(as.character(plot)),
         year = as.numeric(as.character(year)),
         month = as.numeric(as.character(month)),
         day = as.numeric(as.character(day)))

for (i in 1:6) {
  sr_8150[[i]] <- sr_8150[[i]] %>% 
    rename_all(tolower) %>% 
    rename(csr = sr) %>% 
    mutate(block = LETTERS[1]) %>% 
    filter(year %in% c(2012:2013, 2015:2017, 2020), 
           month %in% c(6:9),
           hour %in% c(9:12)) %>% 
    mutate(date = paste(year, month, day, sep = "-")) %>% 
    filter(date %in% sr_8100$date) %>% 
    group_by(block, plot, year, month, day) %>% 
    summarise(csr = mean(as.numeric(csr), na.rm = TRUE),
              .groups = "drop") %>% 
    select(block, plot, year, month, day, csr) %>% 
    na.omit()
}

for (i in 7:12) {
  sr_8150[[i]] <- sr_8150[[i]] %>% 
    rename_all(tolower) %>% 
    rename(csr = sr) %>% 
    mutate(block = LETTERS[6]) %>% 
    filter(year %in% c(2012:2013, 2015:2017, 2020), 
           month %in% c(6:9),
           hour %in% c(9:12)) %>% 
    mutate(date = paste(year, month, day, sep = "-")) %>% 
    filter(date %in% sr_8100$date) %>% 
    group_by(block, plot, year, month, day) %>% 
    summarise(csr = mean(as.numeric(csr), na.rm = TRUE),
              .groups = "drop") %>% 
    select(block, plot, year, month, day, csr) %>% 
    na.omit()
}

sr_8150 <- do.call(rbind, sr_8150) %>% 
  mutate(date = paste(year, month, day, sep = "-")) %>% 
  as.data.frame() %>% 
  mutate(plot = as.numeric(as.character(plot)),
         year = as.numeric(as.character(year)),
         month = as.numeric(as.character(month)),
         day = as.numeric(as.character(day)))

sr <- sr_8100 %>% 
  bind_rows(., sr_8150)


###########################################################
# Generate data -----------------------


# 1. Calculate forage production
# 1.1. Calculate total forage production
forage_production <- spp_tab_clean %>% 
  group_by(year, block, plot) %>% 
  summarise(forage_prod = sum(biomass), 
            .groups = "drop") %>% 
  filter(year %in% year) %>% 
  mutate(exp_period = ifelse(year %in% c(2012, 2013), "Early",
                          ifelse(year == 2020, "Late", "Middle"))) %>% 
  group_by(exp_period, block, plot) %>% 
  summarise(forage_prod = mean(forage_prod),
            .groups = "drop") %>% 
  filter(block != "winter")


# 1.2. Calculate forage production for each plant functional group
forage_production_group <- spp_tab_clean %>% 
  inner_join(., spp_code, by = "species_name") %>% 
  group_by(year, block, plot, functional_group) %>% 
  summarise(group_biomass = sum(biomass),
            .groups = "drop") %>% 
  filter(year %in% year) %>% 
  mutate(exp_period = ifelse(year %in% c(2012, 2013), "Early",
                             ifelse(year == 2020, "Late", "Middle"))) %>% 
  group_by(exp_period, block, plot, functional_group) %>%
  summarise(group_biomass = mean(group_biomass),
            .groups = "drop") %>% 
  filter(block != "winter") %>% 
  pivot_wider(names_from = "functional_group",
              values_from = "group_biomass")


# 2. Calculate forage quality
# 2.1. % productivity accounted by those species with leaf N data
total_biomass <- spp_tab_clean %>% 
  summarise(biomass = sum(biomass)) %>% 
  pull()
leaf_N_species_biomass <- spp_tab_clean %>% 
  filter(species_name %in% leaf_N$species) %>% 
  summarise(biomass = sum(biomass)) %>% 
  pull()
leaf_N_species_biomass/total_biomass*100
# 81.5%

# 2.2 calculate community weighted mean of leaf nitrogen content
spp_com <- spp_tab_clean %>% 
  filter(species_name %in% leaf_N$species) %>% 
  group_by(idd) %>% 
  filter(sum(biomass) > 0) %>% 
  ungroup() %>% 
  group_by(species_name) %>% 
  filter(sum(biomass) > 0) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "species_name",
              values_from = "biomass") %>% 
  select(year, block, plot, levels(factor(leaf_N$species)))

rownames(leaf_N) <- leaf_N$species
leaf_N <- leaf_N[-1]
CWM_leaf_N <- dbFD(leaf_N, spp_com[-c(1:3)],
     w.abun = TRUE, calc.CWM = TRUE,
     calc.FDiv = FALSE, calc.FRic = FALSE)$CWM

forage_quality <- cbind(spp_com[1:3], CWM_leaf_N) %>% 
  filter(year %in% year) %>% 
  mutate(exp_period = ifelse(year %in% c(2012, 2013), "Early",
                             ifelse(year == 2020, "Late", "Middle"))) %>% 
  group_by(exp_period, block, plot) %>%
  summarise(forage_quality = mean(leaf_N),
            .groups = "drop") %>% 
  filter(block != "winter")


# 3. Calculate species richness and beta-diversity (LCBD)
# 3.1. calculate species richness 
species_richness <- spp_tab_clean %>% 
  group_by(year, block, plot, idd) %>% 
  summarise(SR = sum(biomass > 0),
            .groups = "drop")
  

# 3.2. Calculate beta-diversity (LCBD)

spp_com4LCBD <- spp_tab_clean %>% 
  pivot_wider(names_from = "species_name",
               values_from = "biomass") %>% 
  select(Lonicera_rupicola:Carex_pygmaea)
LCBD <- beta.div(spp_com4LCBD,
                 "hellinger")$LCBD
biodiversity <- cbind(species_richness, LCBD) %>%
  filter(year %in% year) %>% 
  mutate(exp_period = ifelse(year %in% c(2012, 2013), "Early",
                             ifelse(year == 2020, "Late", "Middle"))) %>% 
  group_by(exp_period, block, plot) %>%
  summarise(SR = mean(SR),
            LCBD = mean(LCBD),
            .groups = "drop") %>% 
  filter(block != "winter")


# 4. Calculate soil carbon storage

soil_C_early <- 
  soil_BD %>% 
  mutate(layer = forcats::fct_collapse(layer, 
                                       '40-70cm' = c("40-50cm", "50-60cm",
                                                     "60-70cm"))) %>% 
  group_by(block, plot, layer) %>% 
  summarise(BD = mean(BD, na.rm = TRUE),
            .groups = "drop") %>% 
  inner_join(., soil_C_early, by = c("block", "plot", "layer")) %>% 
  mutate(soil_C_stock = ifelse(layer == "40-70cm", 
                               STC*BD*1000*0.3*0.001, 
                               STC*BD*1000*0.1*0.001)) %>% 
  group_by(block, plot) %>% 
  summarise(soil_C_stock = sum(soil_C_stock),
            exp_period = unique(exp_period),
            .groups = "drop")

soil_C_middle <- 
  soil_BD %>% 
  mutate(layer = forcats::fct_collapse(layer, 
                                       '40-70cm' = c("40-50cm", "50-60cm",
                                                     "60-70cm"))) %>% 
  group_by(block, plot, layer) %>% 
  summarise(BD = mean(BD, na.rm = TRUE),
            .groups = "drop") %>% 
  inner_join(., soil_C_middle, by = c("block", "plot", "layer")) %>% 
  mutate(soil_C_stock = ifelse(layer == "40-70cm", 
                               STC*BD*1000*0.3*0.001, 
                               STC*BD*1000*0.1*0.001)) %>% 
  group_by(block, plot) %>% 
  summarise(soil_C_stock = sum(soil_C_stock),
            exp_period = unique(exp_period),
            .groups = "drop")

soil_C_late <- 
  soil_BD %>% 
  inner_join(., soil_C_late, by = c("block", "plot", "layer")) %>% 
  mutate(soil_C_stock = STC*BD*1000*0.1*0.001) %>% 
  group_by(block, plot) %>% 
  summarise(soil_C_stock = sum(soil_C_stock),
            exp_period = unique(exp_period),
            .groups = "drop")

soil_C_stock <- bind_rows(soil_C_early, soil_C_middle) %>% 
  bind_rows(., soil_C_late) %>% 
  mutate(block = tolower(block))


# 5. Calculate soil respiration

soil_respiration <- sr %>% 
  mutate(exp_period = ifelse(year %in% c(2012, 2013), "Early",
                             ifelse(year == 2020, "Late", "Middle"))) %>% 
  group_by(exp_period, block, plot) %>%
  summarise(csr = mean(csr, na.rm = TRUE),
            .groups = "drop") %>% 
  mutate(block = tolower(block))


###########################################################
# combine all the data and save
es_indicator <- merge(forage_production, forage_production_group,
      by = c("exp_period", "block", "plot")) %>% 
  merge(., forage_quality, 
        by = c("exp_period", "block", "plot")) %>% 
  merge(., biodiversity, 
        by = c("exp_period", "block", "plot")) %>% 
  merge(., soil_C_stock, 
        by = c("exp_period", "block", "plot")) %>% 
  merge(., soil_respiration, 
        by = c("exp_period", "block", "plot")) %>% 
  select(-block)

es_indicator <- exp_design %>%
  mutate(plot = as.character(plot)) %>% 
  inner_join(., es_indicator, by = "plot")

saveRDS(es_indicator, "./outputs/ecosystem_service_indicators.rds")


###########################################################
#                   End of the Script                     #
###########################################################