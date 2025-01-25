###########################################################
# Calculating ES weighting
# XJ
# last edited in 2024.04.12
###########################################################

rm(list = ls())

# load library
library(tidyverse)
library(vegan)
library(cowplot)
library(TeachingDemos)
library(ggalluvial)

char2seed("stakeholder")

# load data
questionnaire_data <- read.csv("./data/haibei_questionnaire_data.csv")


# clean data
# number of samples before data cleaning
dim(questionnaire_data)[1]  # 47
questionnaire_data %>% 
  mutate(location = factor(location)) %>% 
  summary()

# since an interviewee may be identified as multiple groups of stakeholders
# we will separate a collapsed column of stakeholder into multiple rows
separate_data <- separate_rows(questionnaire_data, stakeholder) %>% 
  mutate(stakeholder = factor(stakeholder))
summary(separate_data)
# number of samples included
dim(separate_data)[1]  # 60

# coding for stakeholder groups
# A = Environmental protection agency
# B = Biodiversity conservation agency
# C = Climate change mitigation agency
# D =	Policy makers
# E = Pastoralists
# F = Others

# stakeholder groups
summary(factor(separate_data$stakeholder))  # A 14; B 12; C 8; D 7; E 14; other 5

# select A, B, C, E
separate_data <- separate_data %>% 
  filter(stakeholder %in% c("A", "B", "C", "E")) %>% 
  droplevels()

# number of samples after data cleaning
dim(separate_data)[1]  # 48
# number of interviewees
length(unique(separate_data$id))  # 40

separate_data %>% 
  select(id, stakeholder) %>% 
  group_by(id) %>% 
  count() %>% 
  group_by(n) %>% 
  count()

# onsite vs. online # 16 online vs. 24 onsite
separate_data %>% 
  select(id, location) %>% 
  distinct() %>% 
  mutate(location = factor(location)) %>% 
  count(location)

# stakeholder groups
summary(factor(separate_data$stakeholder))  # A 14; B 12; C 8; E 14


###########################################################
# ranking of the ecosystem services in the study area
ranking_data <- separate_data %>% 
  select(stakeholder, q1a:q1f) %>% 
  pivot_longer(cols = q1a:q1f,
               names_to = "variable",
               values_to = "value") %>% 
  filter(variable != "q1d") %>% 
  mutate(stakeholder = factor(stakeholder,
                              levels = c("E", "B", "C", "A"),
                              labels = c("Pastoralist",
                                         "Biodiversity\nconservation agency",
                                         "Climate change\nmitigation agency",
                                         "Environmental\nprotection agency")),
         variable = factor(variable,
                           levels = c("q1a", "q1b", "q1c",
                                      "q1e", "q1f"),
                           labels = c("A place where wildlife lives", 
                                      "A place for\npeople to visit and enjoy", 
                                      "A place to produce food", 
                                      "A place for people to live", 
                                      "A place where\nclimate change is regulated"))) 
overall <- ranking_data %>% 
  group_by(stakeholder, variable) %>% 
  mutate(mu_value = mean(value)) %>% 
  ggplot(aes(x = forcats::fct_reorder(variable, mu_value,
                                      .desc = TRUE), 
             y = value, 
             group = variable)) +
  geom_violin(draw_quantiles = c(0.25, 0.75),
              linetype = "dashed") +
  geom_violin(fill = "transparent", 
              draw_quantiles = 0.5,
              lwd = .5) +
  labs(x = NULL, y = "Ranking") +
  coord_flip() +
  theme_bw() +
  theme(panel.grid = element_blank())
  
ggsave("./outputs/ranking_overall.pdf",
       width = 150, height = 100, units = "mm")


ranking_data %>% 
  ggplot(aes(x = forcats::fct_reorder(variable, value,
                                      .desc = TRUE), 
             y = value, 
             group = variable)) +
  geom_violin(aes(fill = stakeholder)) +
  facet_grid(~ stakeholder) +
  scale_fill_manual(values = c("#00BFC4", "#7CAE00",  
                               "#F8766D", "#bebada")) +
  stat_summary(fun.data = "mean_se", 
               fun.args = list(mult = 1), 
               geom = "pointrange", 
               size = .2,
               color = "black", 
               show.legend = FALSE) + 
  labs(x = NULL, y = "Ranking") +
  coord_flip() +
  theme_bw(base_size = 11.5) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none")

ggsave("./outputs/ranking_by_stakeholder.pdf",
       width = 170, height = 100, units = "mm")


###########################################################
# Stakeholder priorities

# Testing for differences in the structure of stakeholder priorities
# Non-metric multidimensional scaling (NMDS)
stakeholder_weighting <- separate_data %>% 
  select(stakeholder, biodiversity_conservation:co2_emission) %>% 
  mutate(stakeholder = factor(stakeholder,
                              levels = c("E", "B", "C", "A"),
                              labels = c("Pastoralist",
                                         "Biodiversity\nconservation agency",
                                         "Climate change\nmitigation agency",
                                         "Environmental\nprotection agency"))) %>% 
  select(-water_source)  # remove water source

stakeholder_weighting[-1] <- scale(stakeholder_weighting[-1])
dist_matrix <- vegdist(stakeholder_weighting[-c(1, 5)], method = "euclidean")
nmds <- metaMDS(dist_matrix)
nmds_scores <- as.data.frame(scores(nmds))
nmds_scores$stakeholder <- stakeholder_weighting$stakeholder

hull <- nmds_scores %>%
  group_by(stakeholder) %>% 
  slice(chull(NMDS1, NMDS2))

stress_data <- data.frame(x = 3, 
                          y = 3, 
                          stress = glue::glue('Stress = {round(nmds$stress, 2)}'))

overall_nmds <- nmds_scores %>% 
  data.frame() %>% 
  ggplot(aes(NMDS1, NMDS2, color = stakeholder)) +
  geom_vline(xintercept = 0, color = "grey", lwd = .6) +
  geom_hline(yintercept = 0, color = "grey", lwd = .6) +
  geom_polygon(data=hull,
               aes(x=NMDS1,y=NMDS2, 
                   fill = stakeholder, 
                   color = NULL),
               alpha = .5) +
  geom_point(size = 0.8) +
  scale_color_manual(values = c("#00BFC4", "#7CAE00",  
                                "#F8766D", "#bebada"),
                     name = "Stakeholder") +
  scale_fill_manual(values = c("#00BFC4", "#7CAE00",  
                               "#F8766D", "#bebada"),
                    name = "Stakeholder") +
  geom_text(data = stress_data,
            aes(x = x, y = y, label = stress),
            inherit.aes = FALSE) +
  labs(x = "NMDS Axis 1",
       y = "NMDS Axis 2") +
  theme_bw(base_size = 10.5) +
  theme(panel.grid = element_blank(),
        legend.position = "top") +
  guides(color = guide_legend(ncol = 1))

# ggsave("./outputs/NMDS_stakeholder_weighting_remove_water_source.pdf",
#        width = 150, height = 100, units = "mm")


# Analysis of variance using distance matrices
adonis2(dist_matrix ~ stakeholder, 
       data = stakeholder_weighting)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = dist_matrix ~ stakeholder, data = stakeholder_weighting)
#             Df SumOfSqs      R2    F    Pr(>F)   
# stakeholder  3   28.214 0.15007 2.5897  0.013 *
# Residual    44  159.786 0.84993                
# Total       47  188.000 1.00000  


# Multivariate homogeneity of group variances
stakeholder_weighting_betadisper <- betadisper(dist_matrix, 
                                               stakeholder_weighting$stakeholder)
anova(stakeholder_weighting_betadisper)
# Response: Distances
#            Df Sum Sq Mean Sq F value Pr(>F)
# Groups     3  0.429  0.1429  0.1419 0.9343
# Residuals 44 44.322  1.0073  
TukeyHSD(stakeholder_weighting_betadisper)
permutest(stakeholder_weighting_betadisper, 
          pairwise = TRUE, permutations = 999)


###########################################################
# Calculate average stakeholder weightings

# Inverse variance weighting
# https://www.medcalc.org/manual/meta-analysis-generic.php
# In the inverse variance method the weight given to each study is 
# the inverse of the variance of the effect estimate 
# (i.e. one over the square of its standard error). 
# Thus larger studies are given more weight than smaller studies, 
# which have larger standard errors. 
# This choice of weight minimizes the imprecision (uncertainty) of 
# the pooled effect estimate.

stakeholder_weighting <- separate_data %>% 
  select(stakeholder, biodiversity_conservation:co2_emission) %>% 
  select(-water_source) %>% # remove water source
  mutate(stakeholder = factor(stakeholder,
                              levels = c("E", "B", "C", "A"),
                              labels = c("Pastoralist",
                                         "Biodiversity\nconservation agency",
                                         "Climate change\nmitigation agency",
                                         "Environmental\nprotection agency"))) %>% 
  pivot_longer(cols = biodiversity_conservation:co2_emission,
               names_to = "variable", values_to = "value") %>% 
  group_by(stakeholder, variable) %>% 
  summarize(inv_var = (1/sd(value)/sqrt(n()))^2,
            value = mean(value),
            .groups = "drop") %>% 
  mutate(inv_var = ifelse(is.infinite(inv_var), NA, inv_var)) %>% 
  mutate(inv_var = ifelse(is.na(inv_var), mean(inv_var, na.rm = TRUE), inv_var)) %>% 
  mutate(value = value * inv_var) %>% 
  select(-inv_var)

saveRDS(stakeholder_weighting,
        "./outputs/stakeholder_weighting.rds")

###########################################################
# Alluvial plot
# https://r-charts.com/flow/ggalluvial/
alluvial_plot <- stakeholder_weighting %>% 
  # filter(variable != "water_source") %>% 
  # droplevels() %>% 
  mutate(variable = factor(variable,
                           levels = c("forage_production",
                                      "forage_quality",
                                      "biodiversity_conservation",
                                      "carbon_storage", 
                                      "co2_emission"),
                           labels = c("Forage\nproduction",
                                      "Forage\nquality",
                                      "Biodiversity\nconservation",
                                      "Carbon\nstorage", 
                                      "CO2\nmitigation")), # note replace CO2 emission with CO2 mitigation 
         stakeholder = factor(stakeholder,
                              levels = c("Pastoralist",
                                         "Biodiversity\nconservation agency",
                                         "Climate change\nmitigation agency",
                                         "Environmental\nprotection agency"),
                              labels = c("Pastoralist", 
                                         "Biodiversity\nconservation\nagency",
                                         "Climate\nchange\nmitigation\nagency",
                                         "Environmental\nprotection\nagency"))) %>% 
  group_by(stakeholder) %>% 
  mutate(sum_value = sum(value),
         wts = value/sum_value*100) %>% 
  ungroup() %>%
  group_by(variable) %>% 
  mutate(sum_wts = sum(wts),
         rel_wts = wts/sum_wts*100) %>% 
  ggplot(aes(axis1 = stakeholder, axis2 = variable, y = rel_wts)) +
  geom_alluvium(aes(fill = stakeholder)) +
  geom_stratum(width = 1/3) +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("stakeholder", "variable"),
                   expand = c(0.15, 0.05)) +
  scale_fill_manual(values = c("#00BFC4", "#7CAE00",  
                               "#F8766D", "#bebada"),
                    name = "Stakeholder") +
  theme_void(base_size = 10.5) +
  theme(legend.position = "none")

plot_grid(alluvial_plot, overall_nmds,
          ncol = 2,
          rel_widths = c(3/5, 2/5),
          labels = c("(a)", "(b)"))

ggsave("./outputs/stakeholder_weightings.pdf",
       width = 200, height = 120, units = "mm")


###########################################################
##                     End                               ##
###########################################################