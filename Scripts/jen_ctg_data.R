library(tidyverse)

wd <- read_csv("Data/ctg/withdrawal_data.csv") # withdrawal data
trials <- read_csv("Data/ctg/summary_withdrawal_trials.csv") # post 2010 trials, indicating which have withdrawal data
condition <- read_csv("Data/ctg/nct_id_condition.csv") # lookup between SAIL condition and nct_id
cfs2 <- readRDS("Outputs/cfs.rds") %>% filter(term == "(Intercept)") %>% 
  select(c(nct_id, prop_term, sail, atc5)) %>%
  mutate(source = "IPD")

# overall missing data table for single run through clinicaltrials.gov

wd = wd %>% 
  group_by(nct_id) %>% 
  summarise(total_wd = sum(count)) %>%
  ungroup()

# csv created to collect missing data from clinicaltrials.gov
# alldata = alldata %>% mutate(
  # nm_arms = if_else(n_arms_withdrawal != number_of_arms, 1, 0))

# write_csv(alldata, file = "Data/ctg/alldata.csv")

# Manual data check of mismatching number of arms and those with missing data for number of withdrawals
# check column: 0 = assumed to be correct (matching arms, number of withdrawals available)
# check column: 1 = checked and corrected/approved; 2 = no data; 3 = no data but trial published; 

alldata <- read_csv("Data/ctg/alldata.csv") %>%
  filter(checked == 1 | checked == 0)

alldata = alldata %>%
  mutate(prop_term = total_wd / enrollment, 
         atc5 = str_sub(atc_7, 1, 5)) %>%
  rename(sail = index_condition)

#### Heat maps of condition against drug class plotting completion count ####
# remove those with missing completion data (i.e. yoda trials for the moment)
# colour scale represents proportion of early withdrawals: blue = low, red = high
# size represents number of trial participants

pdf("Outputs/heat_map_ctg_count.pdf", width = 18, height = 8)
alldata %>% filter(!is.na(prop_term)) %>% 
  ggplot(aes(x=atc5, y=reorder(sail, desc(sail)), fill = prop_term)) +
  geom_count(aes(colour = prop_term, size = enrollment, alpha = 0.8)) + 
  scale_size_area(max_size = 20) +
  scale_colour_gradient(low = "blue", high = "red") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +
  labs(x="ATC drug class", y="Index condition", colour = "Proportion: early termination", size = "Number of participants") + 
  guides(fill = "none", alpha = "none")
dev.off()

pdf("Outputs/heat_map_ctg_jitter.pdf", width = 18, height = 8)
alldata %>% filter(!is.na(prop_term)) %>% 
  ggplot(aes(x=atc5, y=reorder(sail, desc(sail))), fill = prop_term) +
  geom_jitter(aes(colour = prop_term, size = enrollment, alpha = 0.9)) + 
  scale_size_area(max_size = 20) +
  scale_colour_gradient(low = "blue", high = "red") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +
  labs(x="ATC drug class", y="Index condition", colour = "Proportion: early termination", size = "Number of participants") + 
  guides(alpha = "none")
dev.off()

#### Density plots: IPD vs CTG ####
alldata = alldata %>% mutate(source = "CTG")
cfs_alldata = bind_rows(cfs2, alldata)

pdf("Outputs/density_plot.pdf", width = 12, height = 8)
cfs_alldata %>% 
  ggplot(aes(x=prop_term)) + 
  geom_density(aes(fill = source), alpha = 0.5) + 
  labs(x = "Proportion of early withdrawals", y = "Density") + 
  scale_fill_discrete(name = "Source of trial data", labels = c("CTG", "IPD")) + 
  theme(text = element_text(size = 16))
dev.off()

#### Violin plots: IPD vs CTG ####

alldata = alldata %>% mutate(source = "CTG")
cfs_alldata = bind_rows(cfs2, alldata) %>% distinct(nct_id, .keep_all = TRUE)

pdf("Outputs/violin_plot_all.pdf", width = 20, height = 12)
cfs_alldata %>% filter(!is.na(prop_term)) %>% 
  group_by(nct_id) %>%
  ggplot(aes(x = source, y=prop_term, fill = source)) + 
  geom_violin(scale = "count", draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_jitter(height = 0) +
  theme_bw() +
  theme(text = element_text(size = 20)) +
  labs(x = "Source of trial data", y = "Proportion of early withdrawals") +
  scale_fill_discrete(name = "Source of trial data", labels = c("CTG", "IPD"))
dev.off()

pdf("Outputs/violin_plot_no_jitter.pdf", width = 20, height = 12)
cfs_alldata %>% filter(!is.na(prop_term)) %>% 
  group_by(nct_id) %>%
  ggplot(aes(x = source, y=prop_term, fill = source)) + 
  geom_violin(scale = "count", draw_quantiles = c(0.25, 0.5, 0.75)) +
  theme_bw() +
  theme(text = element_text(size = 20)) +
  labs(x = "Source of trial data", y = "Proportion of early withdrawals") +
  scale_fill_discrete(name = "Source of trial data", labels = c("CTG", "IPD"))
dev.off()
