# Generate research trends in studies of climate change and ecosystem services
# XJ
# 2023/11/18
###########################################################

rm(list = ls())

library(tidyverse)
library(cowplot)


# load data
# data were generated using web of science at Sep 17, 2022
es <- read.table("./data/wos_search/analyze_es.txt",
                 header = TRUE)
es_cc <- read.table("./data/wos_search/analyze_es_cc.txt",
                    header = TRUE)
es_cc_exp <- read.table("./data/wos_search/analyze_es_cc_exp.txt",
                        header = TRUE)
es_cc_exp_long <- read.table("./data/wos_search/analyze_es_cc_exp_long.txt",
                             header = TRUE)
es_cc_exp_long_grass <- read.table("./data/wos_search/analyzees_cc_exp_long_grassland.txt",
                                   header = TRUE)

es_clean <- es %>% 
  bind_rows(., es_cc) %>% 
  bind_rows(., es_cc_exp) %>% 
  bind_rows(., es_cc_exp_long) %>% 
  bind_rows(., es_cc_exp_long_grass) %>% 
  mutate(type = factor(type,
                       levels = c("es", "es_cc", "es_cc_exp",
                                  "es_cc_exp_long", "es_cc_exp_long_grass"),
                       labels = c("Ecosystem services",
                                  "+ Climate change",
                                  "+ Experiment",
                                  "+ Long term",
                                  "+ Grassland")))

# total number of record
es_clean %>% 
  summarise(tot_num = sum(number)) %>% 
  pull(tot_num)
# 62346


###########################################################
# generate plots of research trends
p1 <- es_clean %>% 
  ggplot(aes(year, number, color = type)) +
  geom_line() +
  scale_color_discrete(name = NULL) +
  labs(x = "Year of publication",
       y = "Number of record") +
  theme_classic() +
  theme(legend.position = c(.3, .8))

p2 <- es_clean %>% 
  mutate(rel_value = number/sum(number) * 100) %>% 
  ggplot(aes(year, rel_value, color = type)) +
  geom_line() +
  scale_color_discrete(name = NULL) +
  labs(x = "Year of publication",
       y = "Percentage (%)") +
  theme_classic()

es_rel_clean <- es_clean %>% 
  group_by(type) %>% 
  summarise(tot_value = sum(number)) %>% 
  mutate(tot_rel_value = tot_value/sum(tot_value) * 100) 

p3 <- es_rel_clean %>% 
  ggplot(aes(type, tot_rel_value)) +
  geom_bar(stat = "identity") +
  geom_label(aes(label = round(tot_rel_value, 1))) +
  labs(x = NULL,
       y = "Percentage (%)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, vjust = .5))

plot_grid(p1, p3, nrow = 2,
          # labels = c("(a)", "(b)"),
          align = "v")

ggsave("./outputs/publication_trends.pdf",
       width = 160, height = 160, units = "mm")

###########################################################
#                   End of the Script                     #
###########################################################