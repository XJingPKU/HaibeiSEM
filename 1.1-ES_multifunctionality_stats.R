###########################################################
# Effects of climate change and time on ES and ES-MF
# XJ
# last edited in 2024/04/12
###########################################################

rm(list = ls()) 
# load library
library(tidyverse)
library(cowplot)
library(lmerTest)
library(piecewiseSEM)
library(TeachingDemos)
library(Hmisc)
library(corrplot)
library(RColorBrewer)

char2seed("multifunctionality")

# load data
# df2 <- readRDS("./outputs/es_mf_linear.rds")
df2 <- readRDS("./outputs/es_mf_nonlinear.rds")  # ln-transformed ES indicators
funct_group <- readRDS("./outputs/functional_group_biomass.rds")

df2 <- df2 %>% 
  rename(f.prod = forage_production,
         f.qual = forage_quality,
         b.cons = biodiversity_conservation,
         c.stor = carbon_storage,
         c.flux = co2_emission) %>% 
  mutate(f.prov = 1/2*(f.prod + f.qual),
         c.clim = 1/2*(c.stor + c.flux)) %>% 
  mutate(warm = factor(warm, levels = c("Ambient", "Warming")),
         precip = factor(precip, levels = c("Ambient", "Drought", "Wet")))


###########################################################
# 1. general linear models

# 1.1 calculate explained variance ------------------------
SS <- NULL
for (i in c("f.prod", "f.qual", "b.cons", "c.stor", "c.flux",
            "f.prov", "c.clim",
            "NSH", "EPA", "FH", "BCA", "CCMA")) {
  mod <- aov(df2[, i] ~ block + warm*precip*exp_period, data = df2)
  res <- summary(mod)[[1]]$`Sum Sq`
  SS <- cbind(SS, res)
  # print(model.tables(mod, type = "means", se = TRUE))
}

# variance explained
SS <- data.frame(SS)
names(SS) <- c("f.prod", "f.qual", "b.cons", "c.stor", "c.flux", 
               "f.prov", "c.clim",
               "NSH", "EPA", "FH", "BCA", "CCMA")
rownames(SS) <- c("Block", "Warming", "Precipitation", "Time",
                  "WxP", "WxT", "PxT", "WxPxT", "Resid")
SS <- apply(SS, 2, function(x) {x/sum(x)*100})
SS <- data.frame(SS)
SS$treat <- rownames(SS)
write.csv(SS, "./outputs/variance_explained.csv")


# 1.2 plot explained variance -----------------------------
p.single_ES <- SS %>% 
  select("f.prod", "f.qual", "b.cons", "c.stor", "c.flux", 
         "f.prov", "c.clim", treat) %>% 
  pivot_longer(cols = c(-treat),
               names_to = "variable",
               values_to = "value") %>% 
  data.frame() %>% 
  mutate(treat = factor(treat,
                        levels = c("Warming", "Precipitation", 
                                   "Time",
                                   "WxP", "WxT", "PxT", "WxPxT", "Block", "Resid"),
                        labels = c("Warming (W)", "Precipitation (P)",
                                   "Time (T)", "WxP", "WxT", "PxT", "WxPxT", 
                                   "Block", "Residuals"))) %>% 
  mutate(variable = factor(variable,
                           levels = c("c.clim", "f.prod", "f.qual",
                                      "b.cons", 
                                      "f.prov", "c.stor", "c.flux"),
                           labels = c("Climate change\nregulation",
                                      "Forage production",
                                      "Forage quality",
                                      "Biodiversity\nconservation",
                                      "Forage\nprovision", 
                                      "Carbon storage", 
                                      "CO2 mitigation"))) %>% 
  filter(!treat %in% c("WxPxT", "Block", "Residuals"),
         variable %in% c("Carbon storage", 
                         "CO2 mitigation",
                         "Biodiversity\nconservation",
                         "Forage production",
                         "Forage quality")) %>% 
  ggplot() +
  geom_bar(aes(x = forcats::fct_rev(variable),
               y = value,  
               fill = treat), 
           alpha = 0.9, 
           stat="identity",
           position = position_stack(reverse = TRUE)) +
  labs(x = "",
       y = "Variance explained (%)") +
  ylim(c(0, 95)) +
  scale_fill_manual(values = c("#FC8D62", "#66C2A5", "#8DA0CB", 
                               "#E78AC3",
                               "#A6D854", "#3182bd")) +
  coord_flip() +
  theme_classic(base_size = 12) +
  theme(legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.75, 0.20)) +
  guides(fill = guide_legend(ncol = 2))


p.multi_ES <- SS %>% 
  select("NSH", "EPA", "FH", "BCA", "CCMA", treat) %>% 
  pivot_longer(cols = c(-treat),
               names_to = "variable",
               values_to = "value") %>% 
  data.frame() %>% 
  mutate(treat = factor(treat,
                        levels = c("Warming", "Precipitation", 
                                   "Time",
                                   "WxP", "WxT", "PxT", "WxPxT", "Block", "Resid"),
                        labels = c("Warming (W)", "Precipitation (P)",
                                   "Time (T)", "WxP", "WxT", "PxT", "WxPxT", 
                                   "Block", "Residuals"))) %>% 
  mutate(variable = factor(variable,
                           levels = c("FH", "BCA", "CCMA", "EPA", "NSH"),
                           labels = c("Pastoralist", 
                                      "Biodiversity\nconservation agency", 
                                      "Climate change\nmitigation agency",
                                      "Environmental\nprotection agency",
                                      "Equal weight\nmultifunctionality"))) %>%
  filter(!treat %in% c("WxPxT", "Block", "Residuals")) %>% 
  ggplot() +
  geom_bar(aes(x = forcats::fct_rev(variable),
               y = value,  
               fill = treat), 
           alpha = 0.9, 
           stat="identity",
           position = position_stack(reverse = TRUE)) +
  labs(x = "",
       y = "Variance explained (%)") +
  ylim(c(0, 95)) +
  scale_fill_manual(values = c("#FC8D62", "#66C2A5", "#8DA0CB", 
                               "#E78AC3",
                               "#A6D854", "#3182bd")) +
  coord_flip() +
  theme_classic(base_size = 12) +
  theme(legend.title = element_blank(),
        legend.position = "none")

plot_grid(p.multi_ES, p.single_ES, ncol = 1, align = "v",
          labels = c("(a)", "(b)"))
ggsave("./outputs/es_var_explained.pdf", 
       width = 180, height = 190, units = "mm")


# 1.3 linear mixed-effects models -------------------------

for (i in c("c.stor", "c.flux", "b.cons", "f.prod", "f.qual", 
            "NSH", "CCMA", "BCA", "FH", "EPA")) {
  mod1 <- lmer(df2[, i] ~ warm*precip*exp_period + 
                (1|block/plot), 
               REML = FALSE,
              data = df2)
  mod2 <- update(mod1, .~. -warm:precip:exp_period)
  mod3 <- update(mod2, .~. -precip:exp_period)
  mod4 <- update(mod3, .~. -warm:exp_period)
  mod5 <- update(mod4, .~. -warm:precip)
  mod6 <- update(mod5, .~. -exp_period)
  mod7 <- update(mod6, .~. -precip)
  mod8 <- update(mod7, .~. -warm)
  cat("========================\n", i)
  print(anova(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8))
}


# 1.4 Error bar plot

df2.single <- df2 %>% 
  data.frame() %>% 
  select(exp_period, warm, precip, 
         c.stor, c.flux, b.cons, f.prod, f.qual) %>% 
  mutate(exp_period = factor(exp_period, levels = c("Early", "Middle", "Late")),
         precip = factor(precip, levels = c("Drought", "Ambient", "Wet"))) %>% 
  pivot_longer(cols = c.stor:f.qual,
               names_to = "variable",
               values_to = "value") %>% 
  mutate(variable = factor(variable,
                           levels = c("f.prod", "f.qual", "b.cons", 
                                      "c.stor", "c.flux"),
                           labels = c("Forage production",
                                      "Forage quality",
                                      "Biodiversity\nconservation",
                                      "Carbon storage", 
                                      "CO2 mitigation"))) %>% 
  group_by(exp_period, warm, precip, variable) %>% 
  summarise(mu = mean(value),
            se = sd(value)/sqrt(length(value)),
            .groups = "drop") 


p.single <- df2.single %>% 
  ggplot(aes(exp_period, mu, color = warm, group = warm), alpha = .6) +
  geom_line(lwd = 0.5) +
  geom_point(aes(shape = warm), size = 3, alpha = 0.66) +
  geom_errorbar(aes(ymin = mu - se,
                    ymax = mu + se), 
                # color = "black", 
                alpha = 0.66,
                width = 0.1) +
  facet_grid(variable ~ precip, switch = "y") +
  scale_y_continuous(breaks = c(0.2, 0.6, 1)) +
  expand_limits(y = c(0.2, 1.1)) +
  scale_color_manual(values = c("#0070C0", "red"), name = "") +
  scale_shape_manual(values = c(1, 19), name = "") +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 11.5) +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        panel.grid = element_blank(),
        legend.position = "top")


# egg::tag_facet(p.single)
ggsave("./outputs/es_benefit.pdf", 
       width = 160, height = 185, units = "mm")

# measures of ecosystem multi-serviceability
p.multi <- df2 %>% 
  select(exp_period, warm, precip, NSH, EPA, FH, BCA, CCMA) %>% 
  mutate(exp_period = factor(exp_period, levels = c("Early", "Middle", "Late")),
         precip = factor(precip, levels = c("Drought", "Ambient", "Wet"))) %>% 
  pivot_longer(cols = NSH:CCMA,
               names_to = "variable",
               values_to = "value") %>%  
  mutate(variable = factor(variable,
                           levels = c("FH", "BCA", "CCMA", "EPA", "NSH"),
                           labels = c("Pastoralist",
                                      "Biodiversity\nconservation agency",
                                      "Climate change\nmitigation agency",
                                      "Environmental\nprotection agency",
                                      "Equal weight\nmultifunctionality"))) %>%
  group_by(exp_period, warm, precip, variable) %>% 
  summarise(mu = mean(value),
            se = sd(value)/sqrt(length(value)),
            .groups = "drop") %>% 
  ggplot(aes(exp_period, mu, color = warm, group = warm), alpha = .6) +
  geom_line(lwd = 0.5) +
  geom_point(aes(shape = warm), size = 3, alpha = 0.66) +
  geom_errorbar(aes(ymin = mu - se,
                    ymax = mu + se), 
                # color = "black", 
                alpha = 0.66,
                width = 0.1) +
  facet_grid(variable ~ precip, switch = "y") +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  expand_limits(y = c(0, 1.1)) +
  scale_color_manual(values = c("#0070C0", "red"), name = "") +
  scale_shape_manual(values = c(1, 19), name = "") +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 11.5) +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        panel.grid = element_blank(),
        legend.position = "top")
# egg::tag_facet(p.multi)


ggsave("./outputs/es_mf.pdf", 
       width = 160, height = 185, units = "mm")





###########################################################
##                     End                               ##
###########################################################