# jen_meta_analysis_post.R
library(tidyverse)
library(ggplot2)
library(brms)
library(tidybayes)
library(modelr)
library(RColorBrewer)
library(data.table)

# This will successfully draw a plot with distributions by factor
# .width argument can be changed to display different probability intervals
# dotted reference line OR 1.0, red dashed lines are lower and upper limits from population estimates

cfs <- readRDS("Outputs/cfs.rds")
cfs_como = cfs %>% filter(term == "como_cnt" & model_n == 1 & std.error < 5)
cfs_age = cfs %>% filter(term == "age_std" & model_n == 5 & std.error < 5)
cfs_sex = cfs %>% filter((term == "sex" | term == "sexTRUE") & model_n == 5 & std.error < 5)


# como plots ----

pdf("Outputs/como_sail_plot_with_pooled.pdf")
mod1_brm_como <- readRDS("Scratch_data/condition_model.Rds")
mod1_brm_como <- mod1_brm_como$mdl

a <- mod1_brm_como %>%
  recover_types(cfs_como) %>%
  spread_draws(b_Intercept, r_sail[sail,]) %>%
  mean_qi(sail_mean = exp(b_Intercept + r_sail), .width = c(0.95, 0.8, 0.5)) 
b <- mod0_brm_como %>%
  recover_types(cfs_como) %>%
  spread_draws(b_Intercept) %>%
  mean_qi(sail_mean = exp(b_Intercept), .width = c(0.95, 0.8, 0.5)) %>%
  mutate(sail = "TOTAL ESTIMATE")
ab <- bind_rows(a, b)
ggplot(ab, aes(y=factor(sail, sort(unique(sail), decreasing = TRUE)), x = sail_mean, xmin = .lower, xmax = .upper,
               colour = sail)) +
  geom_pointinterval() +
  scale_x_continuous(limits = c(0.5, 1.75)) +
  scale_colour_manual(values = c("TOTAL ESTIMATE" = "red")) +
  geom_vline(xintercept = 1.0, linetype = "dotted", color = "black") +  
  theme_classic() + 
  labs(x = "OR (50%, 80% and 95% credibility intervals)", y = "Condition") + 
  guides(colour = "none")
dev.off()

pdf("Outputs/como_atc5_plot_with_pooled.pdf")
mod2_brm_como <- readRDS("Scratch_data/drugclass_model.Rds")
mod2_brm_como <- mod2_brm_como$mdl

a <- mod2_brm_como %>%
  recover_types(cfs_como) %>%
  spread_draws(b_Intercept, r_atc5[atc5,]) %>%
  mean_qi(atc5_mean = exp(b_Intercept + r_atc5), .width = c(0.95, 0.8, 0.5)) 
b <- mod0_brm_como %>%
  recover_types(cfs_como) %>%
  spread_draws(b_Intercept) %>%
  mean_qi(atc5_mean = exp(b_Intercept), .width = c(0.95, 0.8, 0.5)) %>%
  mutate(atc5 = "TOTAL ESTIMATE")
ab <- bind_rows(a, b)
ggplot(ab, aes(y=factor(atc5, sort(unique(atc5), decreasing = TRUE)), x = atc5_mean, xmin = .lower, xmax = .upper,
               colour = atc5)) +
  geom_pointinterval() +
  scale_x_continuous(limits = c(0.5, 1.75)) +
  scale_colour_manual(values = c("TOTAL ESTIMATE" = "red")) +
  geom_vline(xintercept = 1.0, linetype = "dotted", color = "black") +  
  theme_classic() + 
  labs(x = "OR (50%, 80% and 95% credibility intervals)", y = "ATC drug class")  + 
  guides(colour = "none")
dev.off()

# age plots ----

pdf("Outputs/age_sail_plot_with_pooled.pdf")
mod1_brm_age <- readRDS("Scratch_data/condition_model_age.Rds")

a <- mod1_brm_age %>%
  recover_types(cfs_age) %>%
  spread_draws(b_Intercept, r_sail[sail,]) %>%
  mean_qi(sail_mean = exp(b_Intercept + r_sail), .width = c(0.95, 0.8, 0.5)) 
b <- mod0_brm_age %>%
  recover_types(cfs_age) %>%
  spread_draws(b_Intercept) %>%
  mean_qi(sail_mean = exp(b_Intercept), .width = c(0.95, 0.8, 0.5)) %>%
  mutate(sail = "TOTAL ESTIMATE")
ab <- bind_rows(a, b)
ggplot(ab, aes(y=factor(sail, sort(unique(sail), decreasing = TRUE)), x = sail_mean, xmin = .lower, xmax = .upper,
               colour = sail)) +
  geom_pointinterval() +
  scale_x_continuous(limits = c(0.5, 1.75)) +
  scale_colour_manual(values = c("TOTAL ESTIMATE" = "red")) +
  geom_vline(xintercept = 1.0, linetype = "dotted", color = "black") +  
  theme_classic() + 
  labs(x = "OR (50%, 80% and 95% credibility intervals)", y = "Condition") + 
  guides(colour = "none")
dev.off()

pdf("Outputs/age_atc5_plot_with_pooled.pdf")
mod2_brm_age <- readRDS("Scratch_data/drugclass_model_age.Rds")

a <- mod2_brm_age %>%
  recover_types(cfs_age) %>%
  spread_draws(b_Intercept, r_atc5[atc5,]) %>%
  mean_qi(atc5_mean = exp(b_Intercept + r_atc5), .width = c(0.95, 0.8, 0.5)) 
b <- mod0_brm_age %>%
  recover_types(cfs_age) %>%
  spread_draws(b_Intercept) %>%
  mean_qi(atc5_mean = exp(b_Intercept), .width = c(0.95, 0.8, 0.5)) %>%
  mutate(atc5 = "TOTAL ESTIMATE")
ab <- bind_rows(a, b)
ggplot(ab, aes(y=factor(atc5, sort(unique(atc5), decreasing = TRUE)), x = atc5_mean, xmin = .lower, xmax = .upper,
               colour = atc5)) +
  geom_pointinterval() +
  scale_x_continuous(limits = c(0.5, 1.75)) +
  scale_colour_manual(values = c("TOTAL ESTIMATE" = "red")) +
  geom_vline(xintercept = 1.0, linetype = "dotted", color = "black") +  
  theme_classic() + 
  labs(x = "OR (50%, 80% and 95% credibility intervals)", y = "ATC drug class")  + 
  guides(colour = "none")
dev.off()

# sex plots ----

pdf("Outputs/sex_sail_plot_with_pooled.pdf")
mod1_brm_sex <- readRDS("Scratch_data/condition_model_sex.Rds")
mod0_brm_sex <- readRDS("Scratch_data/pooled_model_sex.Rds")
cfs_sex <- readRDS("Outputs/cfs_sex.rds")
a <- mod1_brm_sex %>%
  recover_types(cfs_sex) %>%
  spread_draws(b_Intercept, r_sail[sail,]) %>%
  mean_qi(sail_mean = exp(b_Intercept + r_sail), .width = c(0.95, 0.8, 0.5)) 
b <- mod0_brm_sex %>%
  recover_types(cfs_sex) %>%
  spread_draws(b_Intercept) %>%
  mean_qi(sail_mean = exp(b_Intercept), .width = c(0.95, 0.8, 0.5)) %>%
  mutate(sail = "TOTAL ESTIMATE")
ab <- bind_rows(a, b)
ggplot(ab, aes(y=factor(sail, sort(unique(sail), decreasing = TRUE)), x = sail_mean, xmin = .lower, xmax = .upper,
               colour = sail)) +
  geom_pointinterval() +
  scale_x_continuous(limits = c(0.5, 1.75)) +
  scale_colour_manual(values = c("TOTAL ESTIMATE" = "red")) +
  geom_vline(xintercept = 1.0, linetype = "dotted", color = "black") +  
  theme_classic() + 
  labs(x = "OR (50%, 80% and 95% credibility intervals)", y = "Condition") + 
  guides(colour = "none")
dev.off()

pdf("Outputs/sex_atc5_plot_with_pooled.pdf")
mod2_brm_sex <- readRDS("Scratch_data/drugclass_model_sex.Rds")
mod0_brm_sex <- readRDS("Scratch_data/pooled_model_sex.Rds")
cfs_sex <- readRDS("Outputs/cfs_sex.rds")
a <- mod2_brm_sex %>%
  recover_types(cfs_sex) %>%
  spread_draws(b_Intercept, r_atc5[atc5,]) %>%
  mean_qi(atc5_mean = exp(b_Intercept + r_atc5), .width = c(0.95, 0.8, 0.5)) 
b <- mod0_brm_sex %>%
  recover_types(cfs_sex) %>%
  spread_draws(b_Intercept) %>%
  mean_qi(atc5_mean = exp(b_Intercept), .width = c(0.95, 0.8, 0.5)) %>%
  mutate(atc5 = "TOTAL ESTIMATE")
ab <- bind_rows(a, b)
ggplot(ab, aes(y=factor(atc5, sort(unique(atc5), decreasing = TRUE)), x = atc5_mean, xmin = .lower, xmax = .upper,
               colour = atc5)) +
  geom_pointinterval() +
  scale_x_continuous(limits = c(0.5, 1.75)) +
  scale_colour_manual(values = c("TOTAL ESTIMATE" = "red")) +
  geom_vline(xintercept = 1.0, linetype = "dotted", color = "black") +  
  theme_classic() + 
  labs(x = "OR (50%, 80% and 95% credibility intervals)", y = "ATC drug class")  + 
  guides(colour = "none")
dev.off()

# heat maps of condition against drug class plotting completion count ----
# remove those with missing completion data (i.e. yoda trials for the moment)
# colour scale represents proportion of early withdrawals: blue = low, red = high
# size represents number of trial participants

pdf("Outputs/heat_map_ipd_count.pdf", width = 12, height = 8)
cfs <- readRDS("Outputs/cfs.rds")
cfs %>% filter(!is.na(prop_term)) %>% filter(term == "(Intercept)") %>%
  ggplot(aes(x=atc5, y=reorder(sail, desc(sail)), fill = prop_term)) +
  geom_count(aes(colour = prop_term, size = no_participants, alpha = 0.8)) + 
  scale_size_area(max_size = 50) +
  scale_colour_gradient(low = "blue", high = "red") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1), text = element_text(size = 16)) + 
  labs(x="ATC drug class", y="Condition", colour = "Proportion: early termination", size = "Number of participants") + 
  guides(fill = "none", alpha = "none")
dev.off()

pdf("Outputs/heat_map_ipd_jitter.pdf", width = 12, height = 8)
cfs <- readRDS("Outputs/cfs.rds")
cfs %>% filter(!is.na(prop_term)) %>% filter(term == "(Intercept)") %>%
  ggplot(aes(x=atc5, y=reorder(sail, desc(sail))), fill = prop_term) +
  geom_jitter(aes(colour = prop_term, size = no_participants, alpha = 0.9)) + 
  scale_size_area(max_size = 50) +
  scale_colour_gradient(low = "blue", high = "red") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1), text = element_text(size = 16)) +
  labs(x="ATC drug class", y="Condition", colour = "Proportion: early termination", size = "Number of participants") + 
  guides(alpha = "none")
dev.off()

