###########################################################
# Compute ES multifunctionality
# XJ
# last edited in 2024.04.12
###########################################################

rm(list = ls())

# load library
library(tidyverse)
library(vegan)

# load data
es_indicator <- readRDS("./outputs/ecosystem_service_indicators.rds")
stakeholder_weighting <- readRDS("./outputs/stakeholder_weighting.rds")

# clean data
es_indicator_clean <- es_indicator %>% 
  rename_all(tolower) %>% 
  rename(precip = pre,
         forage_qual = forage_quality,
         species_rich = sr,
         soil_resp = csr) %>% 
  mutate(block = factor(block),
         plot = factor(plot),
         exp_period = factor(exp_period,
                       levels = c("Early", "Middle", "Late")),
         warm = factor(warm,
                       levels = c("A", "E"),
                       labels = c("Ambient", "Warming")),
         precip = factor(precip,
                       levels = c("D", "K", "W"),
                       labels = c("Drought", "Ambient", "Wet"))) %>% 
  mutate(treat = factor(paste(warm, precip, sep = "_"),
                        levels = c("Ambient_Ambient", "Warming_Ambient",
                                   "Ambient_Drought", "Ambient_Wet",
                                   "Warming_Drought", "Warming_Wet"),
                        labels = c("Ambient", "Warming", "Drought", 
                                   "Wet", "Warming_Drought", "Warming_Wet"))) %>% 
  select(exp_period, block, plot, treat, warm, precip, 
         forb, grass, legume, sedge,
         forage_prod, forage_qual,
         species_rich, lcbd,
         soil_c_stock, soil_resp) %>% 
  mutate(soil_resp = -1 * soil_resp + max(soil_resp)) %>%  # convert CO2 emission to mitigation
  mutate(log_forb = log(forb),
         log_grass = log(grass),
         log_legume = log(legume),
         log_sedge = log(sedge),
         log_forage_prod = log(forage_prod),
         log_forage_qual = log(forage_qual),
         log_species_rich = log(species_rich),
         log_lcbd = log(lcbd),
         log_soil_c_stock = log(soil_c_stock),
         log_soil_resp = log(soil_resp + 0.001))

# split the data into three datasets
# functional group
funct_group <- es_indicator_clean %>% 
  select(exp_period, block, plot, treat, warm, precip,
         forb, grass, legume, sedge)
funct_group$func_div <- diversity(funct_group[c("forb", "grass", "legume", "sedge")],
                                 index = "shannon")
# saveRDS(funct_group, "./outputs/functional_group_biomass.rds")

# linear ES indicators
es_indicator_linear <- es_indicator_clean %>% 
  select(exp_period, block, plot, treat, warm, precip,
         forage_prod, forage_qual,
         species_rich, lcbd,
         soil_c_stock, soil_resp)

# ln-transformed ES indicators
es_indicator_nonlinear <- es_indicator_clean %>% 
  select(exp_period, block, plot, treat, warm, precip,
         log_forage_prod, log_forage_qual,
         log_species_rich, log_lcbd,
         log_soil_c_stock, log_soil_resp)

es_indicator_linear %>% 
  pivot_longer(cols = forage_prod:soil_resp,
               names_to = "variable",
               values_to = "value") %>% 
  ggplot(aes(exp_period, value)) +
  geom_boxplot() +
  facet_wrap(~ variable, scales = "free")

es_indicator_nonlinear %>% 
  pivot_longer(cols = log_forage_prod:log_soil_resp,
               names_to = "variable",
               values_to = "value") %>% 
  ggplot(aes(exp_period, value)) +
  geom_boxplot() +
  facet_wrap(~ variable, scales = "free")


###########################################################
# Data standardization and ES-MF calculation
###########################################################
# 1. Convert ES indicators into measures of ES supply
rangeStd <- function(x, n.sort) {
  xmin <- mean(sort(x, decreasing = FALSE)[1:n.sort], na.rm = TRUE)
  xmax <- mean(sort(x, decreasing = TRUE)[1:n.sort], na.rm = TRUE)
  x <- (x - xmin) / (xmax - xmin)
  x[which(x > 1)] <- 1
  x[which(x < 0)] <- 0
  return(x)
}

es_supply_linear <- es_indicator_linear %>% 
  mutate(std_forage_prod = rangeStd(forage_prod, 3),  # Standardize ES indicators
         std_forage_qual = rangeStd(forage_qual, 3),
         std_species_rich = rangeStd(species_rich, 3),
         std_lcbd = rangeStd(lcbd, 3),
         std_soil_c_stock = rangeStd(soil_c_stock, 3),
         std_soil_resp = rangeStd(soil_resp, 3)) %>% 
  mutate(forage_production = std_forage_prod,
         forage_quality = std_forage_qual,
         biodiversity_conservation = std_species_rich, 
         carbon_storage = std_soil_c_stock,
         co2_emission = std_soil_resp) %>%  # note CO2 emission should be read as CO2 mitigation
  select(exp_period, block, plot, treat, warm, precip, 
         forage_production, forage_quality,
         biodiversity_conservation,
         carbon_storage, co2_emission)

es_supply_nonlinear <- es_indicator_nonlinear %>% 
  mutate(std_forage_prod = rangeStd(log_forage_prod, 3),  # Standardize ES indicators
         std_forage_qual = rangeStd(log_forage_qual, 3),
         std_species_rich = rangeStd(log_species_rich, 3),
         std_lcbd = rangeStd(log_lcbd, 3),
         std_soil_c_stock = rangeStd(log_soil_c_stock, 3),
         std_soil_resp = rangeStd(log_soil_resp, 3)) %>% 
  mutate(forage_production = std_forage_prod,
         forage_quality = std_forage_qual,
         biodiversity_conservation = std_species_rich, 
         carbon_storage = std_soil_c_stock,
         co2_emission = std_soil_resp) %>% 
  select(exp_period, block, plot, treat, warm, precip, 
         forage_production, forage_quality,
         biodiversity_conservation,
         carbon_storage, co2_emission)

# 2. Convert ES supply into measures of ES-benefit
# 1:1 linear function
es_benefit_linear <- es_supply_linear
es_benefit_nonlinear <- es_supply_nonlinear

# 3. Compute ES-MF
levels(stakeholder_weighting$stakeholder)
stakeholder_weighting_wider <- stakeholder_weighting %>% 
  pivot_wider(names_from = variable,
              values_from = value) %>% 
  select(stakeholder, forage_production,
         forage_quality, 
         biodiversity_conservation,
         carbon_storage, 
         co2_emission) %>% 
  data.frame()
rownames(stakeholder_weighting_wider) <- c("FH", "BCA", "CCMA", "EPA")
stakeholder_weighting_wider <- stakeholder_weighting_wider %>% 
  select(-stakeholder)


es_mf_linear <- es_benefit_linear %>% 
  mutate(NSH = 1/5 * (forage_production + forage_quality + 
                        biodiversity_conservation + 
                        carbon_storage + co2_emission),  # non-stakeholder
         EPA = (forage_production*stakeholder_weighting_wider["EPA", "forage_production"] +
                  forage_quality*stakeholder_weighting_wider["EPA", "forage_quality"] +
                  biodiversity_conservation*stakeholder_weighting_wider["EPA", "biodiversity_conservation"] +
                  carbon_storage*stakeholder_weighting_wider["EPA", "carbon_storage"] +
                  co2_emission*stakeholder_weighting_wider["EPA", "co2_emission"])/
           sum(stakeholder_weighting_wider["EPA", ]),  # Environmental protection
         BCA = (forage_production*stakeholder_weighting_wider["BCA", "forage_production"] +
                  forage_quality*stakeholder_weighting_wider["BCA", "forage_quality"] +
                  biodiversity_conservation*stakeholder_weighting_wider["BCA", "biodiversity_conservation"] +
                  carbon_storage*stakeholder_weighting_wider["BCA", "carbon_storage"] +
                  co2_emission*stakeholder_weighting_wider["BCA", "co2_emission"])/
           sum(stakeholder_weighting_wider["BCA", ]),  # Biodiversity conservation
         CCMA = (forage_production*stakeholder_weighting_wider["CCMA", "forage_production"] +
                   forage_quality*stakeholder_weighting_wider["CCMA", "forage_quality"] +
                   biodiversity_conservation*stakeholder_weighting_wider["CCMA", "biodiversity_conservation"] +
                   carbon_storage*stakeholder_weighting_wider["CCMA", "carbon_storage"] +
                   co2_emission*stakeholder_weighting_wider["CCMA", "co2_emission"])/
           sum(stakeholder_weighting_wider["CCMA", ]), # Climate change mitigation
         FH = (forage_production*stakeholder_weighting_wider["FH", "forage_production"] +
                 forage_quality*stakeholder_weighting_wider["FH", "forage_quality"] +
                 biodiversity_conservation*stakeholder_weighting_wider["FH", "biodiversity_conservation"] +
                 carbon_storage*stakeholder_weighting_wider["FH", "carbon_storage"] +
                 co2_emission*stakeholder_weighting_wider["FH", "co2_emission"])/
           sum(stakeholder_weighting_wider["FH", ]))  # Farmer&herder
         

es_mf_nonlinear <- es_benefit_nonlinear %>% 
  mutate(NSH = 1/5 * (forage_production + forage_quality + 
                        biodiversity_conservation + 
                        carbon_storage + co2_emission),  # non-stakeholder
         EPA = (forage_production*stakeholder_weighting_wider["EPA", "forage_production"] +
                  forage_quality*stakeholder_weighting_wider["EPA", "forage_quality"] +
                  biodiversity_conservation*stakeholder_weighting_wider["EPA", "biodiversity_conservation"] +
                  carbon_storage*stakeholder_weighting_wider["EPA", "carbon_storage"] +
                  co2_emission*stakeholder_weighting_wider["EPA", "co2_emission"])/
           sum(stakeholder_weighting_wider["EPA", ]),  # Environmental protection
         BCA = (forage_production*stakeholder_weighting_wider["BCA", "forage_production"] +
                  forage_quality*stakeholder_weighting_wider["BCA", "forage_quality"] +
                  biodiversity_conservation*stakeholder_weighting_wider["BCA", "biodiversity_conservation"] +
                  carbon_storage*stakeholder_weighting_wider["BCA", "carbon_storage"] +
                  co2_emission*stakeholder_weighting_wider["BCA", "co2_emission"])/
           sum(stakeholder_weighting_wider["BCA", ]),  # Biodiversity conservation
         CCMA = (forage_production*stakeholder_weighting_wider["CCMA", "forage_production"] +
                   forage_quality*stakeholder_weighting_wider["CCMA", "forage_quality"] +
                   biodiversity_conservation*stakeholder_weighting_wider["CCMA", "biodiversity_conservation"] +
                   carbon_storage*stakeholder_weighting_wider["CCMA", "carbon_storage"] +
                   co2_emission*stakeholder_weighting_wider["CCMA", "co2_emission"])/
           sum(stakeholder_weighting_wider["CCMA", ]), # Climate change mitigation
         FH = (forage_production*stakeholder_weighting_wider["FH", "forage_production"] +
                 forage_quality*stakeholder_weighting_wider["FH", "forage_quality"] +
                 biodiversity_conservation*stakeholder_weighting_wider["FH", "biodiversity_conservation"] +
                 carbon_storage*stakeholder_weighting_wider["FH", "carbon_storage"] +
                 co2_emission*stakeholder_weighting_wider["FH", "co2_emission"])/
           sum(stakeholder_weighting_wider["FH", ]))  # Farmer&herder

# save data
saveRDS(funct_group, "./outputs/functional_group_biomass.rds")
saveRDS(es_mf_linear, "./outputs/es_mf_linear.rds")
saveRDS(es_mf_nonlinear, "./outputs/es_mf_nonlinear.rds")


###########################################################
##                     End                               ##
###########################################################